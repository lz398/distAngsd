#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <chrono>
#include<htslib/bgzf.h>
#include<htslib/kstring.h>
#include "shared.h"
#include <pthread.h>
#include "GLtest.h"
using namespace std;
using namespace Eigen;

// Test Code which should be removed later on;
void tic()
{
    t0 = Clock::now();
}
void toc()
{
    Clock::time_point t1 = Clock::now();
    milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
    std::cout <<"Elapsed time is "<< ms.count() << " milliseconds\n";
}

/* Optimization function */
// Brent's method, Golden section search for the minimum of f(x)
#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double brent(double ax, double bx, double cx, double (*f)(double), double tol,  double *xmin)
{
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0;
    
    a=(ax < cx ? ax : cx); // a=min{ax,cx}
    b=(ax > cx ? ax : cx); // b=max{ax,cx}
    x=w=v=bx; // All stay on the right end of the interval; Keeps x values of previous 2 steps
    fw=fv=fx=(*f)(x); // All equal to the function value at the right end of the interval. Keeps funtion values corresponding to x of previous 2 steps
    for (int iter=1;iter<=ITMAX;iter++) {
        xm=0.5*(a+b); // midpoint of (a,b)=(ax,cx)
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
        } // When the interval is small enough, regard as convergence
        if (fabs(e) > tol1) {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            else {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        } else {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(u);
        if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }
    printf("Too many iterations in brent");
    *xmin=x;
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI

void SeedSetup(){
    long int seed = -1;
    //  ++argv;
    //  while(*argv){
    //    if(strcasecmp(*argv,"-seed")==0)
    //      seed = atol(*(++argv));
    //    /*
    //      else if(strcasecmp(*arg,"-flag"))
    //      flag =0;
    //    */
    //    else
    //      fprintf(stderr,"\t-> Unknown argument: %s\n",*argv);
    //    ++argv;
    //  }
    seed = time(NULL) ^getpid();
    fprintf(stderr,"seed: %ld\n",seed);
    //cout << "Random seed is "<< seed<<"\n";
    //seed it
    //first version rand family
    srand48(seed);
    //second version drand48 family
    //srand48(seed);
}

void SetSeed(int seed)
{
    z_rndu = 170*(seed%178) + 137;
}

// gamma distribution
double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947,-86.50532032942, 24.01409824083,-1.231739572450, 0.1208650973866e-2,-0.5395239385e-5};
    
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -=(x+0.5)*log(tmp);
    ser=1.00000000019015;
    for (j=0; j<5; j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// random problem
// Poisson distribution
double Poisson(double xm) {
    //srand48(int(time));
    double gammln(double xx);
    static double sq,alxm,g,oldm=(-1.0);
    double em, t, y;
    
    if (xm < 12.0) {
        if (xm != oldm) {
            oldm=xm;
            g=exp(-xm);
        }
        em=-1;
        t=1.0;
        do {
            ++em;
            t *=drand48();
        } while (t>g);
    }
    else {
        if (xm!=oldm) {
            oldm=xm;
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                y=tan(M_PI*drand48());
                em=sq*y+xm;
            } while (em< 0.0);
            em=floor(em);
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
        } while (drand48()>t);
    }
    return em;
}


double uniform()
{
    /*
     U(0,1): AS 183: Appl. Stat. 31:188-190
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190
     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.
     
     Suggested to me by Ziheng Yang who also provided me with
     the source code used here.  I use it because it is both fast and portable.
     */
    static int x_rndu=11, y_rndu=23;
    double r;
    
    x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
    y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
    z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
    if (x_rndu<0) x_rndu+=30269;
    if (y_rndu<0) y_rndu+=30307;
    if (z_rndu<0) z_rndu+=30323;
    r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
    return (r-(int)r);
}

/*Different substitution models*/
//hard coded transition probabilities under JC69 model.  '*par' is not used as JC models uses no parameters beyond t.
double pijtJC(double t, double *par, int from, int to)
{
    if (from==-1 || to==-1) return 1.0; /*GAP  - can be deleted if sequences have no gaps*/
    if (t>150.0) return 0.25;
    else if (from==to) return 0.25+0.75*exp(-4.0*t/3.0);
    else return 0.25-0.25*exp(-4.0*t/3.0);
}

//hard coded transition probabilities under HKY85 model
//This code has not been tested
double pijtHKY(double t, double *par, int from, int to)
{      /*ACGT*/
    double A, PIj, kappa;
    
    if (from==-1 || to==-1) return 1.0; /*GAP  - can be deleted if sequences have no gaps*/
    if (t>150.0) return par[to];
    kappa=par[4];
    if  (to==0 || to==2) PIj=par[0]+par[2];
    else PIj=par[1]+par[3];
    A=1.0+PIj*(kappa-1.0);
    if (from==to)
        return par[to]+par[to]*exp(-t)*(1.0/PIj-1.0)+exp(-A)*((PIj-par[to])/PIj);
    else if (from+to==2 || from+to==4)
        return par[to]+par[to]*(1.0/PIj-1.0)*exp(-t)-(par[to]/PIj)*exp(-t*A);
    else return par[to]*(1.0-exp(-t));
}

//hard coded transition probabilities under F84 model
//this code has not been tested
double pijtF84(double t, double *par, int from, int to)
{      /*ACGT*/
    double A, PIj, kappa;
    
    if (from==-1 || to==-1) return 1.0; /*GAP  - can be deleted if sequences have no gaps*/
    if (t>150.0) return par[to];
    kappa=par[4];
    if  (to==0 || to==2) PIj=par[0]+par[2];
    else PIj=par[1]+par[3];
    A=kappa+1.0;
    if (from==to)
        return par[to]+par[to]*exp(-t)*(1.0/PIj-1.0)+exp(-A)*((PIj-par[to])/PIj);
    else if (from+to==2 || from+to==4)
        return par[to]+par[to]*(1.0/PIj-1.0)*exp(-t)-(par[to]/PIj)*exp(-t*A);
    else return par[to]*(1.0-exp(-t));
}

void diagonalizeGTR( double *par)
/*parameters are
 0: A<->C
 1: A<->G
 2: A<->T
 3: C<->G
 4: C<->T
 (G <->T is defined as 1.0 to avoid identifiability issues)
 5: piA
 6: piC
 7:piG
 (pA equals 1-piA-piC-piG)
 */
{
    double stand;
    
    MatrixXd A(4,4);
    A(0,0)= 0.0; A(0,1) =  par[0];   A(0,2) = par[1]; A(0,3) = par[2];
    A(1,0) = par[0]; A(1,1) =  0.0;  A(1,2) = par[3]; A(1,3) = par[4];
    A(2,0) = par[1]; A(2,1) =  par[3];   A(2,2) = 0.0; A(2,3) = 1.0;
    A(3,0) = par[2]; A(3,1) =  par[4];   A(3,2) = 1.0; A(3,3) = 0.0;
    
    pi(0)=par[5]; pi(1)=par[6]; pi(2)=par[7]; pi(3)=1.0-pi(0)-pi(1)-pi(2);
    if (pi(3)<0) {printf("Pi error in maketransitionmatrix"); exit(-1);}
    
    stand=0.0;
    for (int i=0; i<4;i++){
        for (int j=0; j<4; j++){
            if (i!=j) {
                A(i,j)=A(i,j)*pi(j);
                A(i,i)=A(i,i)-A(i,j);
            }
        }
        stand=stand-pi(i)*A(i,i);
    }
    for (int i=0; i<4;i++){
        for (int j=0; j<4; j++){
            A(i,j)=A(i,j)/stand;
        }
    }
    
    /*printf("par: ");
     for (i=0; i<8;i++)
     printf("%lf ",par[i]);
     printf("\n");
     
     for (i=0; i<4;i++){
     for (j=0; j<4; j++)
     printf("%lf ",A[i][j]);
     printf("\n");
     }*/
    
    EigenSolver<MatrixXd> es(A);
    VectorXcd RVAL(4);
    VectorXd  RIVAL(4);
    MatrixXcd RVEC(4,4);
    MatrixXd  RIVEC(4,4);
    RVAL=es.eigenvalues();
    RRVAL=RVAL.real();
    RIVAL=RVAL.imag();
    RVEC=es.eigenvectors();
    RRVEC=RVEC.real();
    RIVEC=RVEC.imag();
    
    if (abs(RRVAL.maxCoeff())>1e-8 || RIVAL.norm()>1e-8){
        cout<<RRVAL;
        printf("Transitions matrix did not converge or contained non-real values!\n");
        exit(-1);
    }
    LRVEC=RRVEC;
    
    if (RRVEC.determinant()==0){
        printf("Could not invert matrix!\nResults may not be reliable!\n");
        exit(-1);
    }else{
        RRVEC=RRVEC.inverse();
    }
    //if (eigen(1, A[0], 4, RRVAL, RIVAL, RRVEC[0], RIVEC[0], workspace) != 0){
    //    printf("Transitions matrix did not converge or contained non-real values!\n");
    //    exit(-1);
    //    }
    //for (i=0; i<4; i++)
    //    for (j=0; j<4; j++)
    //        LRVEC[i][j] = RRVEC[i][j];
    //if (matinv(RRVEC[0],4, 4, workspace) != 0){
    //    printf("Could not invert matrix!\nResults may not be reliable!\n");
    //    exit(-1);
    //    }
    newmat=0;
}

void gettransitionprobmatGTR(double t)
{
    VectorXd EXPOS(4);
    int i, j, k;
    
    for (k=0; k<4; k++)
    EXPOS(k) = exp(t*RRVAL(k));
    for (i=0; i<4; i++){
        for (j=0; j<4; j++){
            PMAT(i,j) = 0.0;
            for (k=0; k<4; k++)
            PMAT(i,j) =  PMAT(i,j) + RRVEC(k,j)*LRVEC(i,k)*EXPOS(k);
        }
    }
    newt=0;
    /*printf("PMAT:\n");
     for (i=0; i<4;i++){
     for (j=0; j<4; j++)
     printf("%lf ",PMAT[i][j]);
     printf("\n");
     }*/
    
    
}

double pijtGTR(double t, double *par,  int from, int to)
{
    
    if (t>150.0) return pi(to);
    if (newmat==1) {diagonalizeGTR(par); gettransitionprobmatGTR(t);}
    else if (newt==1) gettransitionprobmatGTR(t);
    if (from==-1 || to==-1) return 1.0; /*GAP  - can be deleted if sequences have no gaps*/
    else return PMAT(from,to);
}

//double like(int numsites, int **data, double (*pijt)(double t, double *par, int from, int to), double *locpi, double t, double *par)
//{
//    int i, j;
//    double L=0;
//
//    newmat=newt=1;/*these variables set to 1 specifies that the transition rate matrix needs to be diagonalized and that the transition prob matrix needs to be recalculated*/
//    /*note that this code needs to be moved out of the function if JC or HKY functionality is implemented*/
//    for (i=0; i<numsites; i++){
//        j=data[0][i];
//        L += log(locpi[j]*pijt(t,par,j,data[1][i]));
//    }
//    return L;
//}

//Might be put into io.h
//int readseqnuc()
//{
//    FILE *infile;
//    int i, j, numsites;
//    double PIstandard;
//    char c;
//
//    if (NULL==(infile=fopen("nuc.txt","r"))){
//        puts ("Cannot find the infile named nuc.txt");
//        exit(-1);}
//    fscanf(infile,"%i",&numsites);
//    DATA = (int**) malloc(2*(sizeof(int *)));
//    for (i=0; i<2; i++)
//    DATA[i]=(int*) malloc(numsites*(sizeof(int)));
//    while ((c=(fgetc(infile)))!='\n');
//    i=pi(0)=pi(1)=pi(2)=pi(3)=0;
//    do{
//        j=0;
//        while ((c=fgetc(infile))!='\n')
//        {
//            if ((c!=' ')&&(c!='\t'))
//            {
//                if (j>=numsites){printf("ERROR READING DATA!!\n");  scanf("%i",&i); exit(-1);}
//                if (c=='a' || c=='A')
//                {
//                    DATA[i][j] = 0;
//                    pi(0)++;
//                }
//                else if (c=='c'||c=='C')
//                {
//                    DATA[i][j] = 1;
//                    pi(1)++;
//                }
//                else if (c=='g'||c=='G')
//                {
//                    DATA[i][j] = 2;
//                    pi(2)++;;
//                }
//                else if (c=='t' || c=='u' || c=='T' || c== 'U')
//                {
//                    DATA[i][j] = 3;
//                    pi(3)++;
//                }
//                else if (c=='n' || c=='-' || c=='N')
//                    DATA[i][j] = -1;
//                else{ printf("\nBAD BASE in sequence %i base %i: %c",i+1,j+1,c);  scanf("%i",&i); exit(-1);}
//                j++;
//            }
//        }
//        i++;
//    } while (i<2);
//    PIstandard = 0.0;
//    for (i=0; i<4; i++)
//    PIstandard+=pi(i);
//    for (i=0; i<4; i++)
//    pi(i)=pi(i)/PIstandard;
//    fclose(infile);
//    return numsites;
//}

/* General simulation*/
void simnucleotides(int Genotype[2], double simmat[4][4])
{
    double u;
    int i, j;
    
    u = uniform();
    i=j=0;
    while (u>simmat[j][i])
    {
        if (i==3){
            i=0;
            j++;
        }
        else{
            i++;
        }
        if (j>3) {printf("Numerical error 1 in simulation algorithm"); exit(-1);}
    }
    Genotype[0]=i;
    Genotype[1]=j;
}

void makeSIMMAT(double simmat[4][4])
{
    for (int j=0; j< 4; j++){
        if (j==0) simmat[0][0]=pi(0)*PMAT(0,0);
        else simmat[j][0]=pi(j)*PMAT(j,0)+simmat[j-1][3];
        for (int i=1; i<4; i++)
        simmat[j][i]=simmat[j][i-1]+pi(j)*PMAT(j,i);
    }
    if (simmat[3][3]<0.9999999 || simmat[3][3]>1.0000001) {
        printf("Numerical error 2 in simulation algorithm"); exit(-1);
    }else{
        simmat[3][3]=1.0;
    }
}

void makepolyMAT(double simmat[4][4])
{
    for (int j=0; j< 4; j++){
        simmat[j][0]=PMAT(j,0);
        for (int i=1; i<4; i++)
        simmat[j][i]=simmat[j][i-1]+PMAT(j,i);
        if (simmat[j][3]<0.9999999 || simmat[j][3]>1.0000001){
            printf("Numerical error 3 in simulation algorithm"); exit(-1);
        }else{
            simmat[j][3]=1.0;  /*not sure if this costly protection against numerical errors is necessary*/
        }
    }/*
      printf("PMAT:\n");
      for (j=0; j< 4; j++){
      for (i=0; i<4; i++)
      printf("%lf ",PMAT[j][i]);
      printf("\n");
      }
      printf("simmat:\n");
      for (j=0; j< 4; j++){
      for (i=0; i<4; i++)
      printf("%lf ",simmat[j][i]);
      printf("\n");
      }*/
}

void simpoly(int nuc, double SIMMAT[4][4], int genotype[2])
{
    double u;
    for (int i=0;i<2;i++){
        u=uniform();
        int j=0;
        while (u>SIMMAT[nuc][j]){
            j++;
            if (j>3) {
                printf("Numerical error 4 in simulation algorithm"); exit(-1);
            }
        }
        genotype[i]=j;
    }
}

/*we hardcode this to avoid having to do any arithmetic in the conversion and to make sure homozygous types are checked first*/
/*this function is optimized to be fast and not to be pretty*/
/*
 AA=00: 0
 AC=01: 1
 AG=02: 2
 AT=03: 3
 CC=11: 4
 CG=12: 5
 CT=13: 6
 GG=22: 7
 GT=23: 8
 TT=33: 9
 */
int findgenotypeindex(int i, int j)
{
    int k;
    
    if (i>j){
        k=i;
        i=j;
        j=k;
    }
    if (j<4){
        if (i==0){
            if (j==0) return 0;
            else if (j==1) return 1;
            else if (j==2) return 2;
            else if (j==3) return 3;
            else {printf("error in genotype index conversion table"); exit(-1);}
        }
        else if (i==1){
            if (j==1) return 4;
            else if (j==2) return 5;
            else if (j==3) return 6;
            else {printf("error in genotype index conversion table"); exit(-1);}
        }
        else if (i==2){
            if (j==2) return 7;
            else if (j==3) return 8;
            else {printf("error in genotype index conversion table"); exit(-1);}
        }
        else if (i==3 && j==3) return 9;
        else {printf("error in genotype index conversion table"); exit(-1);}
    }else if(j==4){
        if (i == 0){
            return 10;
        } else if (i==1){
            return 11;
        } else if (i==2){
            return 12;
        } else if (i==3){
            return 13;
        } else if (i==4){
            return 14;
        } else {printf("error in genotype index conversion table"); exit(-1);}
    }else if (j==5){
        if (i<=5){
            if (i == 0){
                return 10;
            } else if (i==1){
                return 11;
            } else if (i==2){
                return 12;
            } else if (i==3){
                return 13;
            } else if (i==4){
                return 14;
            } else if (i==5){
                return 14;
            }
        } else {printf("error in genotype index conversion table"); exit(-1);}
    }else if(j==6){
        if (i<=6){
            return 15;
        } else {printf("error in genotype index conversion table"); exit(-1);}
    }
    else{printf("error in genotype index conversion table");}
    exit(-1);
}


//int findgenotypeindex(int i, int j)
//{
//    int k;
//
//    if (i>j){
//        k=i;
//        i=j;
//        j=k;
//    }
//    if (j<4){
//        if (i==0){
//            if (j==0) return 0;
//            else if (j==1) return 1;
//            else if (j==2) return 2;
//            else if (j==3) return 3;
//            else {printf("error in genotype index conversion table"); exit(-1);}
//        }
//        else if (i==1){
//            if (j==1) return 4;
//            else if (j==2) return 5;
//            else if (j==3) return 6;
//            else {printf("error in genotype index conversion table"); exit(-1);}
//        }
//        else if (i==2){
//            if (j==2) return 7;
//            else if (j==3) return 8;
//            else {printf("error in genotype index conversion table"); exit(-1);}
//        }
//        else if (i==3 && j==3) return 9;
//        else {printf("error in genotype index conversion table"); exit(-1);}
//    }else if(j==4){
//        if (i == 0){
//            return 10;
//        } else if (i==1){
//            return 11;
//        } else if (i==2){
//            return 12;
//        } else if (i==3){
//            return 13;
//        } else if (i==4){
//            return 14;
//        } else {printf("error in genotype index conversion table"); exit(-1);}
//    }else if (j==5){
//        if (i<=5){
//            if (i == 0){
//                return 10;
//            } else if (i==1){
//                return 11;
//            } else if (i==2){
//                return 12;
//            } else if (i==3){
//                return 13;
//            } else if (i==4){
//                return 15;
//            } else if (i==6){
//                return 14;
//            }
//        } else {printf("error in genotype index conversion table"); exit(-1);}
//    }else if(j==6){
//        if (i<=6){
//            return 16; // We will treat them as no information for now. But may come back to here later. _Lei
//        } else {printf("error in genotype index conversion table"); exit(-1);}
//    }
//    else{printf("error in genotype index conversion table"); exit(-1);}
//
//}
//int findgenotypeindex(int i, int j)
//{
//    int k;
//
//    if (i>j){
//        k=i;
//        i=j;
//        j=k;
//    }
//    if (i==0){
//        if (j==0) return 0;
//        else if (j==1) return 1;
//        else if (j==2) return 2;
//        else if (j==3) return 3;
//        else {printf("error in genotype index conversion table"); exit(-1);}
//    }
//    else if (i==1){
//        if (j==1) return 4;
//        else if (j==2) return 5;
//        else if (j==3) return 6;
//        else {printf("error in genotype index conversion table"); exit(-1);}
//    }
//    else if (i==2){
//        if (j==2) return 7;
//        else if (j==3) return 8;
//        else {printf("error in genotype index conversion table"); exit(-1);}
//    }
//    else if (i==3 && j==3) return 9;
//    else {printf("error in genotype index conversion table"); exit(-1);}
//}

void findgenotypes_from_index(int inde, int genotype[2])
{
    switch(inde) {
        case 0: genotype[0]=0; genotype[1]=0; break;
        case 1: genotype[0]=0; genotype[1]=1; break;
        case 2: genotype[0]=0; genotype[1]=2; break;
        case 3: genotype[0]=0; genotype[1]=3; break;
        case 4: genotype[0]=1; genotype[1]=1; break;
        case 5: genotype[0]=1; genotype[1]=2; break;
        case 6: genotype[0]=1; genotype[1]=3; break;
        case 7: genotype[0]=2; genotype[1]=2; break;
        case 8: genotype[0]=2; genotype[1]=3; break;
        case 9: genotype[0]=3; genotype[1]=3; break;
        default: {printf("error 2 in genotype index conversion table (case=%i)",inde); exit(-1);}
    }
}

//Not used
// This function generates the emission probability/switching probability from the true genotype to the observed genotype. But as the observations are short reads instead of genotypes, and the error rates could vary among different rates even at one base, we will temperly not this function.
/*we assume that all errors are symmetric and independent so that
 switching one base in the genotype has probability e and switching two has probability e^2*/
void makeGLswitchmatrix(double e, double switchmatrix[10][10], double simswitchmatrix[10][10])
{
    for (int i=0; i<4; i++){
        for (int j=i; j<4; j++){
            for (int k=0; k<4; k++){
                for (int v=k; v<4; v++){
                    int m1 = findgenotypeindex(i, j);
                    int m2 = findgenotypeindex(v, k);
                    if (i==j){/*originally homozygous*/
                        if (v==k){
                            if (i==k) switchmatrix[m1][m2]=(1.0-e)*(1.0-e); /*remained in same state*/
                            else switchmatrix[m1][m2]=e*e/9.0; /*switched to another homozygous state*/
                        }
                        else if ((v==i||k==i)||(v==j||k==j)) switchmatrix[m1][m2]=2.0*(1-e)*e/3.0; /*switched to a heterozygous state with 1 nucleotide identical*/
                        else switchmatrix[m1][m2]=2.0*e*e/9.0; /*switched to a heterozygous state with 0 nucleotides identical*/
                    }
                    else {/*originally heterozygous*/
                        if  (i==k && j==v) switchmatrix[m1][m2]=(1.0-e)*(1.0-e) + e*e/9.0; /*stayed the same*/
                        else if (v==k) {
                            if  (i==k || j==k) switchmatrix[m1][m2] = (1-e)*e/3.0; /*switched to homozygous state with one nucleotide identical*/
                            else switchmatrix[m1][m2] = e*e/9.0;   /*switched to homozygous state with zero nucleotides identical*/
                        }
                        else if ((v==i||k==i)||(v==j||k==j)) switchmatrix[m1][m2]=(1-e)*e/3.0 + e*e/9.0; /*switched to heterozygous state with one identical*/
                        else  switchmatrix[m1][m2]= 2.0*e*e/9.0;  /*switched to heterozygous state with none identical*/
                    }
                }
            }
        }
    }
    
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            if (j==0) simswitchmatrix[i][j]=switchmatrix[i][j];
            else simswitchmatrix[i][j]=switchmatrix[i][j]+simswitchmatrix[i][j-1];
        }
    }
    
    /*double sum;
     for (j=0; j<10; j++){
     sum=0.0;
     for (i=0; i<10; i++){ if (i==j) switchmatrix[i][j] = switchmatrix[i][j]*2000; //CHANGE11111
     sum += switchmatrix[i][j];}
     for (i=0; i<10; i++)
     switchmatrix[i][j] = switchmatrix[i][j]/sum;
     }
     for (i=0; i<10; i++){
     for (j=0; j<10; j++)
     printf("%lf ",switchmatrix[i][j]);
     printf("\n");
     }*/
}

int simSEQs_per_read(int genotypes[2], double e)
{
    int out;
    int in = genotypes[(int)floor(2.0*uniform())];
    double u = uniform();
    if (u < e) {
        out = (in+1+(int)floor(3.0*uniform()))%4;
    }else{
        out = in;
    }
    return out;
}

void simSEQs_reads(int **SEQDATA, int genotypes[2], size_t site, int species, double e, double RD)
{
    int ReadDepth = ReaddepthGenerator(RD, 0);
    int nuccount[4]={0,0,0,0};
    int out;
    for(int i=0;i<ReadDepth;i++){
        //cout<<simGLs_per_read(GLDATA, genotypes, site, species, e)<<"\n";
        //double e1 = -e*log(1-uniform());
        out = simSEQs_per_read(genotypes, e);
        nuccount[out] = nuccount[out]+1;
    }
    double prob[4];
    prob[0] = (double)nuccount[0]/(double)ReadDepth;
    int j = 0;
    double u = uniform();
    while(u > prob[j]){
        j = j+1;
        prob[j] = prob[j-1]+(double)nuccount[j]/(double)ReadDepth;
    }
    SEQDATA[site][species]=j;
}



/* 1.2. Joint Genotype distribution method: standard and tree structure*/
//Simulation: Generate a read of a true genotype with calling error rate e, and update the genotype likehood.
int simGLs_per_read(double **GLDATA, int genotypes[2], size_t site, int species, double e)
{
    int out;
    int in = genotypes[(int)(2*uniform())];
    double u = uniform();
    if (u < e) {
        out = (in+1+(int)(3*uniform()))%4;
    }else{
        out = in;
    }
    //cout << "read is " << out <<"\n";
    double prob;
    int gene[2];
    for (int k=0; k<10; k++){
        findgenotypes_from_index(k, gene);
        if (out == gene[0]){
            prob = (1.0 - e)/2.0;
        }else{
            prob = e/6.0;
        }
        if (out == gene[1]){
            prob = prob + (1.0 - e)/2.0;
        }else{
            prob = prob + e/6.0;
        }
        if (prob <=0){
            GLDATA[site][10*species+k] = GLDATA[site][10*species+k] - 1000000.0;
        }else{
            GLDATA[site][10*species+k] = GLDATA[site][10*species+k] + log(prob);
        }
    }
    return out;
}

/* Simulation: Generate Read depth for each site*/
int ReaddepthGenerator(double RD, int par){
    switch(par){
        case 1:{
            int ReadDepth = Poisson(RD);
            while (ReadDepth ==0){
                ReadDepth = Poisson(RD);
            }
            break;
        }
        case 2:{
            return floor(RD);
            break;
        }
        default:{
            return Poisson(RD);
        }
    }
}

// Simulation: Gather different reads and finalize the genotype likelihood for one site.
/* In the future should be able to switch from Poisson distributed read depth to constant coverage,
 and deal with 0 read depth.*/
void simGLs_reads(double **GLDATA, int genotypes[2], size_t site, int species, double e, double RD)
{
    for(int k=0;k<10;k++){
        GLDATA[site][10*species+k] = 0;
    }
    //cout<<"Site "<<site<<" Genotype is "<<genotypes[0]<<genotypes[1]<<"\n";
    int ReadDepth = ReaddepthGenerator(RD, 0);
    //cout << "ReadDepth is " << ReadDepth << "\n";
    for(int i=0;i<ReadDepth;i++){
        //cout<<simGLs_per_read(GLDATA, genotypes, site, species, e)<<"\n";
        //double e1 = -e*log(1-uniform());
        //        double e1 = -e*log(1-uniform());
        //        if (e1 > 0.5) {e1=0.5;}
        //        simGLs_per_read(GLDATA, genotypes, site, species, e1);
        simGLs_per_read(GLDATA, genotypes, site, species, e);
    }
    if (ReadDepth > 0){
        for(int k=0;k<10;k++){
            GLDATA[site][10*species+k] = exp(GLDATA[site][10*species+k]);
        }
    }
}

//Simulation: Simulate two individuals' genome and reads in two species.
/*tdiv is the divergence time. t1 and t2 are the average coalescence times within species*/
void simulateGLsTwoSpecies(double RD, size_t numsites, double errorate, double tdiv, double t1, double t2, double **GLDATA,  double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (size_t i=0; i<numsites; i++)
    ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    int ** TrueGenome = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (size_t i=0; i<numsites; i++)
    TrueGenome[i]=(int *) malloc(4*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simnucleotides(ancDATA[i], SIMMAT);
        // printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    cout<<"Ancient Data was derived!\n";
    
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    //makeGLswitchmatrix(errorate, switchmatrix, simswitchmatrix);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][0], SIMMAT, genotypes);
        //       cout<<genotypes[0]<<"\n";
        TrueGenome[i][0] = genotypes[0];
        TrueGenome[i][1] = genotypes[1];
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        //    simGLs(GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix);
        simGLs_reads(GLDATA, genotypes, i, 0, errorate, RD);
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][1], SIMMAT, genotypes);
        TrueGenome[i][2] = genotypes[0];
        TrueGenome[i][3] = genotypes[1];
        //     simGLs(GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix);
        simGLs_reads(GLDATA, genotypes, i, 1, errorate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    //cout<<TrueGenome[100][0]<<"\n";
//    for (int i=0; i<1000; i++){
//        for (int j=0; j<4; j++){
//            cout << TrueGenome[i][j] << "\t";
//        }
//        cout << "\n";
//    }
//    int d1 = 0;
//    int d2 = 0;
//    int d3 = 0;
//    for (int s = 0; s<numsites; s++){
//        for (int i = 0; i < 2; i++){
//            for (int j = 2; j < 4; j++){
//                if (TrueGenome[s][i]!=TrueGenome[s][j]){
//                    d3 = d3 + 1;
//                }
//            }
//        }
//        if (TrueGenome[s][0]!=TrueGenome[s][1]){
//            d1 = d1 + 1;
//        }
//        if (TrueGenome[s][2]!=TrueGenome[s][3]){
//            d2 = d2 + 1;
//        }
//    }
//    double l1 = (double)d3*0.25;
//    double l2 = ((double)d1+(double)d2)*0.25;
//    cout << (l1-l2)/l1 <<"\n";
//    cout << (t1+t2)/2*l1/l2<<"\n";
//    cout << (t1+t2)/2/(1+log(1-(l1-l2)/l1))<<"\n";
    //cout << -(t1+t2)/2*0.75*log(1-(l1-l2)/l1*(double)4/(double)3)<<"\n";
    //cout << -log(1-(l1-l2)/l1) <<"\n";
    /*then we free memory*/
    for (size_t i=0; i<numsites; i++)
    free(ancDATA[i]);
    free(ancDATA);
}

//Filter: All sites considered should have data for both individual.
size_t FilterSites(double **GLDATA, int *SDATA, size_t numsites){
    size_t eff_numsites = 0;
    for (size_t s=0; s<numsites; s++){
        int k = 0;
        int l1 = 0;
        int l2 = 0;
        while (l1==0 && k<10){
            if (GLDATA[s][k]!=0){
                l1 = 1;
            }
            k = k+1;
        }
        k = 0;
        while (l2==0 && k<10){
            if (GLDATA[s][10+k]!=0){
                l2 = 1;
            }
            k = k+1;
        }
        SDATA[s] = l1*l2;
        if (SDATA[s]>0){
            eff_numsites = eff_numsites + 1;
        }
    }
    cout << "Effective number of sites is "<<eff_numsites<<"\n";
    return eff_numsites;
}

//Inference: The following two functions are used to measure the convergence of joint-distribution of genotypes.
void differr2D(double* m1, double* m2, double* mdiff, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            mdiff[i*n+j] = m1[i*n+j] - m2[i*n+j];
        }
    }
}

double err2D(double* m, int n){
    double e = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            e = e + pow(m[i*n+j],2);
        }
    }
    return e;
}

//Inference: One EM step for joint genotype distribution
//numsites need to be size_t!
void EMStepfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int *SDATA, size_t numsites){
    double p[10][10], sum;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            ESFS2[i][j]=0.0;
        }
    }
    for (size_t s=0; s<numsites; s++){
        if (SDATA[s]>0){
            sum=0;
            for (int i=0; i<10; i++){
                for (int j=0; j<10; j++){
                    p[i][j] = twoDSFS[i][j]*GLDATA[s][i]*GLDATA[s][10+j];
                    sum += p[i][j];
                }
            }
            for (int i=0; i<10; i++)
            for (int j=0; j<10; j++)
            ESFS2[i][j] += p[i][j]/sum;
        }
    }
}


 void *athread(void *aptr){
     EMjob* ptr = (EMjob*) aptr;
     size_t index = ptr->index;
     double** gls = ptr->gls;
     int *ss = ptr->ss;
     size_t start = ptr -> start;
     size_t len = ptr->len;
     int rowlen = ptr->rowlen;
 //    double sum = 0;
 //    for (;at<=start+len-1;at++){
 //        double *p_at = gls[at];
 //        for (int k = 0; k < rowlen; k++){
 //            sum += *(p_at+k);
 //        }
 //    }
     EMStepfor2DSFS(ptr->two2DSFS, ptr->segESFS2, &gls[start], &ss[start], len);
 //    ptr->sum = sum;
 //    fprintf(stderr,"[%s] thread %lu: Position from %lu to %lu: Finished\n",__FUNCTION__,index,start,start+len-1);
     return NULL;
 }
 
void EMStepfor2DSFS_threading_initial(size_t numsites,int rowL,int nthreads,vector<EMjob> &jobvec, double **GLDATA, int *SDATA){
     size_t n = (floor)((double)numsites/(double)nthreads);
     EMjob tmp;
     tmp.index = 0;
     tmp.gls = GLDATA;
     tmp.ss = SDATA;
     tmp.start = 0;
     tmp.len = n;
     tmp.rowlen = rowL;
     for (int j=0;j<10;j++){
         for (int k=0;k<10;k++){
             tmp.two2DSFS[j][k] = 0.01;
         }
     }
     jobvec.push_back(tmp);
     for (size_t k=1;k<nthreads-1;k++){
         tmp.index = k;
         //tmp.gls =GLDATA;
         //tmp.p_start = testvec.back().p_start+n;
         tmp.start = jobvec.back().start+n;
         tmp.len = n;
         tmp.rowlen = rowL;
         jobvec.push_back(tmp);
     }
     tmp.index = nthreads-1;
     //tmp.p_start = testvec.back().p_start+n;
     //tmp.gls =GLDATA;
     tmp.start = jobvec.back().start+n;
     tmp.len = numsites-(nthreads-1)*n;
     tmp.rowlen = rowL;
     jobvec.push_back(tmp);
     //
     //cout<<testvec.size();
     fprintf(stderr,"[%s] beginning\n",__FUNCTION__);
 }
 
 void EMStepfor2DSFS_threading(pthread_t* mythd, int nthreads,vector<EMjob> &jobvec,double **GLDATA,double two2DSFS[10][10], size_t numsites, size_t eff_numsites){
     size_t i,j,k;
     for(i=0;i<nthreads;i++){
         int rc = pthread_create(&mythd[i],NULL,athread,&jobvec[i]);
         if(rc)
             fprintf(stderr,"problem starting threadid: %lu\n",i);
     }
     for(i=0;i<nthreads;i++){
         if(pthread_join(mythd[i], NULL)) {
             fprintf(stderr, "Error joining thread\n");
         }
     }
     for (i=0; i<10; i++){
         for (j=0; j<10; j++){
             two2DSFS[i][j]=0.0;
         }
     }
     for (j=0;j<10;j++){
         for (k=0;k<10;k++){
             for(i=0;i<nthreads;i++){
                 two2DSFS[j][k] += jobvec[i].segESFS2[j][k];
             }
             two2DSFS[j][k] = two2DSFS[j][k]/(double)eff_numsites;
             for(i=0;i<nthreads;i++){
                 jobvec[i].two2DSFS[j][k] = two2DSFS[j][k];
             }
         }
     }
 }

 

// Inference: Accelerate EM step for joint genotype distribution
int EMAccelfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int *SDATA, size_t numsites){
    double EM1[10][10], EM2[10][10];
    double stepMin =1;
    double stepMax0 = 1;
    double mstep=4;
    static double stepMax=stepMax0;
    EMStepfor2DSFS(twoDSFS, EM1, GLDATA, SDATA, numsites);
    double differr1[10][10];
    differr2D(&EM1[0][0], &twoDSFS[0][0], &differr1[0][0],10);
    double err1 = err2D(&differr1[0][0],10);
    if(sqrt(err1)<tole){
        for(int i=0;i<10;i++){
            for(int j=0;j<10;j++){
                ESFS2[i][j] = EM1[i][j];
            }
        }
        return 1;
    }
    
    EMStepfor2DSFS(EM1, EM2, GLDATA, SDATA, numsites);
    double differr2[10][10];
    differr2D(&EM2[0][0], &EM1[0][0], &differr2[0][0], 10);
    double err2 = err2D(&differr2[0][0], 10);
    if(sqrt(err2)<tole){
        for(int i=0;i<10;i++){
            for(int j=0;j<10;j++){
                ESFS2[i][j] = EM2[i][j];
            }
        }
        return 1;
    }
    
    double differr3[10][10];
    differr2D(&differr2[0][0], &differr1[0][0], &differr3[0][0], 10);
    double err3 = err2D(&differr3[0][0],10);
    double alpha = sqrt(err1/err3);
    alpha = max(stepMin,min(stepMax,alpha));
    for(int i=0;i<10;i++){
        for(int j=0;j<10;j++){
            ESFS2[i][j] = twoDSFS[i][j]+2*alpha*differr1[i][j]+pow(alpha,2)*differr3[i][j];
        }
    }
    
    int outofparspace =0;
    for(int i=0;i<10;i++){
        for(int j=0;j<10;j++){
            if(ESFS2[i][j]<0){
                outofparspace++;
                // break;
            }
        }
    }
    if(outofparspace){
        for(int i=0;i<10;i++){
            for(int j=0;j<10;j++){
                ESFS2[i][j] = EM2[i][j];
            }
        }
    }
    
    if (fabs(alpha - 1) > 0.01){
        //        niter++;
        double ESFS2_tmp[10][10];
        EMStepfor2DSFS(ESFS2,ESFS2_tmp, GLDATA, SDATA, numsites);
        for(int i=0;i<10;i++){
            for(int j=0;j<10;j++){
                swap(ESFS2[i][j],ESFS2_tmp[i][j]);
            }
        }
    }
    //
    if ((alpha - stepMax) > -0.001) {
        stepMax = mstep*stepMax;
    }
    return 1;
}

//Inference: Main EM algorithm for joint genotype distribution
double estimate2DSFS_EM(double twoDSFS[10][10], double **GLDATA, int *SDATA, size_t eff_numsites, size_t numsites)
{
    double ESFS2[10][10], ptemp, d;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            twoDSFS[i][j]=0.01;
        }
    }
    //tic();
    do {
        EMStepfor2DSFS(twoDSFS, ESFS2, GLDATA, SDATA, numsites);
        //EMAccelfor2DSFS(twoDSFS, ESFS2, GLDATA, numsites);
        d=0.0;
        for (int i=0; i<10; i++){
            for (int j=0; j<10; j++){
                ptemp = ESFS2[i][j]/(double)eff_numsites;
                d += pow(ptemp-twoDSFS[i][j],2);
                twoDSFS[i][j]=ptemp;
            }
        }
        //printf("d: %lf\n",d);
    }while (sqrt(d)> tole);
    //toc();
    return d;
}


double estimate2DSFS_EM_threading(double twoDSFS[10][10], double **GLDATA, int* SDATA, size_t numsites, size_t eff_numsites ,int nthreads, int rowL)
 {
     double ptemp[10][10], d;
     for (int i=0; i<10; i++){
         for (int j=0; j<10; j++){
             twoDSFS[i][j]=0.01;
         }
     }
     
     vector<EMjob> jobvec;
     EMStepfor2DSFS_threading_initial(numsites,rowL,nthreads,jobvec,GLDATA,SDATA);
     pthread_t *mythd = new pthread_t[nthreads];
     cout << jobvec.size()<<"\n";
     
     //tic();
     do {
         EMStepfor2DSFS_threading(mythd,nthreads,jobvec,GLDATA,ptemp,numsites,eff_numsites);
         //EMStepfor2DSFS(twoDSFS, ESFS2, GLDATA, numsites);
         //EMAccelfor2DSFS(twoDSFS, ESFS2, GLDATA, numsites);
         d=0.0;
         for (int i=0; i<10; i++){
             for (int j=0; j<10; j++){
                 d += pow(ptemp[i][j]-twoDSFS[i][j],2);
                 twoDSFS[i][j]=ptemp[i][j];
             }
         }
         //printf("d: %lf\n",d);
     }while (sqrt(d)> tole);
     
     delete [] mythd;
     //toc();
     return d;
 }
 

double seqlike(int observedbase, int truebase, double errorrate)
{
    if (observedbase==truebase) return 1.0 - errorrate;
    else return errorrate/3.0;
}

//Inference: total likelihood based on the inferred joint genotype distribution.
double likeGLwithtwoDSFS(double twoDSFS[10][10], double t, double par[8])
{
    int g1[2], g2[2];
    double like, totlike=0.0;
    //cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            like =0.0;
            for (int k=0; k<2; k++){
                for (int v=0; v<2; v++)
                like += log(pi(g1[k])*PMAT(g1[k],g2[v]));
            }
            totlike += like*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

//Inference: total likelihood based on the inferred joint genotype distribution.
//Add in the tree structure
double likeGLwithtwoDSFS_m(double twoDSFS[10][10], double t, double par[8], double t1, double t2)
{
    int g1[2], g2[2];
    double like, totlike=0.0;
    Matrix4d PMAT1, PMAT2;
    diagonalizeGTR(par);
    gettransitionprobmatGTR(t1);
    PMAT1=PMAT;
    gettransitionprobmatGTR(t2);
    PMAT2=PMAT;
    
    gettransitionprobmatGTR(t-t1-t2);
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            like =0.0;
            for (int k=0; k<4; k++){
                for (int v=0; v<4; v++){
                    if(g1[0]==g1[1] && g2[0]==g2[1]){
                        like +=pi(k)*pow(PMAT1(k,g1[0]),2)*PMAT(k,v)*pow(PMAT2(v,g2[0]),2);
                    }else if(g1[0]==g1[1]){
                        like +=pi(k)*pow(PMAT1(k,g1[0]),2)*PMAT(k,v)*2*PMAT2(v,g2[0])*PMAT2(v,g2[1]);
                    }else if(g2[0]==g2[1]){
                        like +=pi(k)*2*PMAT1(k,g1[0])*PMAT1(k,g1[1])*PMAT(k,v)*pow(PMAT2(v,g2[0]),2);
                    }else{
                        like +=pi(k)*4*PMAT1(k,g1[0])*PMAT1(k,g1[1])*PMAT(k,v)*PMAT2(v,g2[0])*PMAT2(v,g2[1]);
                    }
                }
            }
            totlike += log(like)*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT(double t)
{
    return likeGLwithtwoDSFS(GLOBtwoDSFS, t,  GLOBpar);
}

//Inference: Calculate the likelihood for divergence t, joint genotype distribution + tree structure
double likelihoodforT_m(double t)
{
    return likeGLwithtwoDSFS_m(GLOBtwoDSFS, t,  GLOBpar, globt1, globt2);
}

//Inference: Estimation of divergence time t based on joint genotype distribution
double estimateT(double twoDSFS[10][10], double *t, double parameters[])
{
    /*OK this is stupid*/
    for(int i=0;i<10;i++){
        for (int j=0; j<10; j++){
            GLOBtwoDSFS[i][j]=twoDSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

/*Inference: Estimation of divergence time t based on joint genotype distribution + tree structure.*/
double estimateT_m(double twoDSFS[10][10], double *t, double parameters[], double t1, double t2)
{
    /*OK this is stupid*/
    for(int i=0;i<10;i++){
        for (int j=0; j<10; j++){
            GLOBtwoDSFS[i][j]=twoDSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    globt1=t1;
    globt2=t2;
    
    double MLV = brent(0.0000001, 1.0, 10.0, likelihoodforT_m, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution*/
double testtwoDSFS(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double par[9], const char* glfname, int isthreading, int dobinary, int r)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
//    pi(0)=0.25;
//    pi(1)=0.25;
//    pi(2)=0.25;
//    pi(3)=0.25;
//    for (int i=0; i<5; i++)
//    parameters[i]=1.0 /*+ (double)i*/;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
    
//    pi(0)=0.2184;
//    pi(1)=0.2606;
//    pi(2)=0.3265;
//    pi(3)=0.1946;
//    parameters[0]=2.0431;
//    parameters[1]=0.0821;
//    parameters[2]=0.0000;
//    parameters[3]=0.0670;
//    parameters[4]=0.0000;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
//
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    if (glfname!=NULL){
        string str = glfname+to_string(r);
        gls_writer_double(str.c_str(), dobinary, numsites, GLDATA);
    }
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            cout << GLDATA[i][j] << "\t";
    //        }
    //        cout << "\n";
    //    }
    //
    
    //Filter the effective numsites
    size_t eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    if (isthreading==1){
        estimate2DSFS_EM_threading(twoDSFS, GLDATA, SDATA, numsites, eff_numsites, 25, 2*10);
    }else{
        estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    }
    
    
    
    //Estimate T
    estimateT(twoDSFS, &t, parameters);
    cout<<"Estimated t = "<<t<<"\n";
    
    
    for (size_t i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);
    
    free(SDATA);
    return t;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + tree structure*/
double testtwoDSFS_m(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (size_t i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
//    pi(0)=0.2184;
//    pi(1)=0.2606;
//    pi(2)=0.3265;
//    pi(3)=0.1946;
//    parameters[0]=2.0431;
//    parameters[1]=0.0821;
//    parameters[2]=0.0000;
//    parameters[3]=0.0670;
//    parameters[4]=0.0000;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //Filter the effective numsites
    size_t eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //Estimate T
    estimateT_m(twoDSFS, &t, parameters,t1,t2);
    cout<<"Estimated t = "<<t<<"\n";
    
    for (size_t i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);
    free(SDATA);
    return t;
}


/* 3. A Random Read method */
/*Simulation: Generate a read of a true genotype with calling error rate e, and based on the chosen read update likelihood.*/
int simSEQs_per_read_v1(double **SEQDATA, int genotypes[2], size_t site, int species, double e, int flag)
{
    int out;
    int in = genotypes[(int)floor(2.0*uniform())];
    
    double u = uniform();
    if (u < e) {
        out = (in+1+(int)floor(3.0*uniform()))%4;
        //cout<<"in: "<<in<<" "<<"out: "<<out<<"\n";
    }else{
        out = in;
    }
    //cout<<"species "<<species<<" site "<<site<<" out :"<<out<<"\n";
    if (flag == 1){
        double a = log(1.0-e);
        double b = log(e)-log(3.0);
        for(int k=0;k<4;k++){
            if (k == out){
                SEQDATA[site][4*species+k] = SEQDATA[site][4*species+k] + a;
            }else{
                SEQDATA[site][4*species+k] = SEQDATA[site][4*species+k] + b;
            }
        }
    }
    return out;
}

/*Simulation: Gather different reads and finalize the chosen read likelihood for one site. In the future should be able to switch from Poisson distributed read depth to constant coverage, and deal with 0 read depth.*/
void simSEQs_reads_v1(double **SEQDATA, int genotypes[2], size_t site, int species, double e, double RD)
{
    int ReadDepth = ReaddepthGenerator(RD, 0);
    int nuccount[4]={0,0,0,0};
    int out;
    for(int k=0;k<4;k++){
        SEQDATA[site][4*species+k] = 0;
    }
    int j=(int)floor((double)ReadDepth*uniform());
    for(int i=0;i<ReadDepth;i++){
        //cout<<simGLs_per_read(GLDATA, genotypes, site, species, e)<<"\n";
        if (i!=j){
            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e,0);
        }else{
            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e,1);
        }
        //        double e1 = -e*log(1-uniform());
        //        if (e1 > 0.5) {e1=0.5;}
        //        if (i!=j){
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,0);
        //        }else{
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,1);
        //        }
        nuccount[out] = nuccount[out]+1;
    }
    if (ReadDepth>0){
        for(int k=0;k<4;k++){
            SEQDATA[site][4*species+k] = exp(SEQDATA[site][4*species+k]);
            //cout<<SEQDATA[site][4*species+k]<<"\t";
        }
    }
    //cout<<"\n";
}


/*Simulation: Simulate two individuals' genome and reads in two species.
 tdiv is the divergence time. t1 and t2 are the average coalescence times within species, based on the chosen read*/
void simulateGLsTwoSpeciesSEQ_v1(double RD, size_t numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (size_t i=0; i<numsites; i++)
    ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simnucleotides(ancDATA[i], SIMMAT);
        //printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    cout<<"Ancient Data was derived!\n";
    
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    // makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][0], SIMMAT, genotypes);
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        // simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_v1(SEQDATA, genotypes, i, 0, errorrate, RD);
        
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][1], SIMMAT, genotypes);
        //simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_v1(SEQDATA, genotypes, i, 1, errorrate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    /*then we free memory*/
    for (size_t i=0; i<numsites; i++)
    free(ancDATA[i]);
    free(ancDATA);
}

//Filter: All sites considered should have data for both individual.
size_t FilterSitesSEQ(double **SEQDATA, int *SEQ_SDATA, size_t numsites){
    size_t eff_numsites = 0;
    for (size_t s=0; s<numsites; s++){
        int k = 0;
        int l1 = 0;
        int l2 = 0;
        while (l1==0 && k<4){
            if (SEQDATA[s][k]!=0){
                l1 = 1;
            }
            k = k+1;
        }
        k = 0;
        while (l2==0 && k<4){
            if (SEQDATA[s][4+k]!=0){
                l2 = 1;
            }
            k = k+1;
        }
        SEQ_SDATA[s] = l1*l2;
        if (SEQ_SDATA[s]>0){
            eff_numsites = eff_numsites + 1;
        }
    }
    cout << "Effective number of sites is "<<eff_numsites<<"\n";
    return eff_numsites;
}

//Inference: total likelihood based on a random chosen/consensus(future) read per site
double SEQlikelihood_v1(double **SEQDATA, int *SEQ_SDATA, double t, double *par, size_t numsites)
{
    double like, likesite;
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    like=0.0;
    for (size_t s=0; s<numsites; s++){
        if (SEQ_SDATA[s] > 0){
            likesite=0.0;
            for (int i=0; i<4; i++){
                for (int j=0; j<4; j++){
                    likesite+=pi(i)*PMAT(i,j)*SEQDATA[s][i]*SEQDATA[s][4+j];
                    //    printf("%i%i: %lf %lf %lf\n",i,j,pi[i]*PMAT[i][j],seqlike(SEQDATA[s][0],i, errorrate),seqlike(SEQDATA[s][1],j,errorrate));
                }}
            like += log(likesite);
        }
    }
    return -like;
}

/*Inference: Calculate the likelihood based on a random chosen/consensus(future) read per site*/
double likelihoodforTSEQ_v1(double t)
{
    return SEQlikelihood_v1(SEQDATA, SEQ_SDATA, t,  GLOBpar, globnumsites);
}

/*Simulation+Inference: Simulation and estimation of divergence time t based on a randomly chosen read per site*/
double testsimSEQDATA_v1(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate)
{
    double t, MLV, parameters[8];
    
    
    //SetSeed(6);
    SetSeed(rand()%30000+1);
    for (int i=0; i<4; i++)
    pi(i)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
//    pi(0)=0.2184;
//    pi(1)=0.2606;
//    pi(2)=0.3265;
//    pi(3)=0.1946;
//    parameters[0]=2.0431;
//    parameters[1]=0.0821;
//    parameters[2]=0.0000;
//    parameters[3]=0.0670;
//    parameters[4]=0.0000;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
    
    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    SEQDATA = (double **) malloc(numsites * sizeof(double *));
    for (size_t i = 0; i < numsites; i++)
    SEQDATA[i] = (double *) malloc(8 * sizeof(double));
    
    SEQ_SDATA = (int *) malloc(numsites * sizeof(int));
    
    simulateGLsTwoSpeciesSEQ_v1(RD, numsites, errorrate,  tdiv,  t1,  t2, SEQDATA, pijtGTR, parameters);
    
    size_t eff_numsites = FilterSitesSEQ(SEQDATA, SEQ_SDATA, numsites);
    
    for (int i=0; i<8; i++)
    GLOBpar[i]=parameters[i];
    globnumsites=numsites;
    globerror=errorrate;
    MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ_v1, 0.000001, &t);
    cout<<"Estimated t = "<<t<<"\n";
    
    for (size_t i = 0; i < numsites; i++)
    free(SEQDATA[i]);
    free(SEQDATA);
    free(SEQ_SDATA);
    return t;
}

/* 4. Joint nucleotide distribution method */
/*Simulation: Generate a read of a true genotype with calling error rate e*/
int simSEQs_per_read_v2(int genotypes[2], size_t site, int species, double e)
{
    int out;
    int in = genotypes[(int)floor(2.0*uniform())];
    
    double u = uniform();
    if (u < e) {
        out = (in+1+(int)floor(3.0*uniform()))%4;
        //cout<<"in: "<<in<<" "<<"out: "<<out<<"\n";
    }else{
        out = in;
    }
    return out;
}

/*Simulation: Gather different reads and calculate nucleotide likelihood per read for one site. In the future should be able to switch from Poisson distributed read depth to constant coverage, and deal with 0 read depth.*/
void simSEQs_reads_v2(vector<vector<double4> > &P, int genotypes[2], size_t site, int species, double e, double RD)
{
    int ReadDepth = ReaddepthGenerator(RD, 0);
    int nuccount[4]={0,0,0,0};
    int out;
    vector<double4> pstr;
    for(int i=0;i<ReadDepth;i++){
        //cout<<simGLs_per_read(GLDATA, genotypes, site, species, e)<<"\n";
        out=simSEQs_per_read_v2(genotypes, site, species, e);
        double4 str;
        for (int j=0;j<4;j++){
            if(out == j){
                str.vec[j] = 1.0-e;
            }else{
                str.vec[j] = e/3.0;
            }
        }
        pstr.push_back(str);
        //        double e1 = -e*log(1-uniform());
        //        if (e1 > 0.5) {e1=0.5;}
        //      out=simSEQs_per_read_v2(genotypes, site, species, e1);
        nuccount[out] = nuccount[out]+1;
    }
    P.push_back(pstr);
    //cout<<"\n";
}

/*Simulation: Simulate two individuals' genome and reads in two species.
 tdiv is the divergence time. t1 and t2 are the average coalescence times within species, based on the joint selected nucleotide distribution*/
void simulateGLsTwoSpeciesSEQ_v2(double RD, size_t numsites, double errorrate, double tdiv, double t1, double t2, vector<vector<double4> >&P0, vector<vector<double4> >&P1, double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (size_t i=0; i<numsites; i++)
    ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simnucleotides(ancDATA[i], SIMMAT);
        //printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    cout<<"Ancient Data was derived!\n";
    
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    // makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][0], SIMMAT, genotypes);
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        // simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_v2(P0, genotypes, i, 0, errorrate, RD);
        
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        simpoly(ancDATA[i][1], SIMMAT, genotypes);
        //simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_v2(P1, genotypes, i, 1, errorrate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    /*then we free memory*/
    for (size_t i=0; i<numsites; i++)
    free(ancDATA[i]);
    free(ancDATA);
}

/*Check the effective number of sites*/
size_t CheckSites(vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites){
    size_t eff_numsites = 0;
    for (int s=0;s<numsites;s++){
        if (P0[s].size()>0 && P1[s].size()>0){
            eff_numsites += 1;
        }
    }
    cout << "Effective number of sites is "<<eff_numsites<<"\n";
    return eff_numsites;
}

/*Inference: total likelihood based on the inferred joint selected nucleotide distribution.*/
double likeSEQwithtwoDSFS(double SEQ2DSFS[4][4], double t, double par[8])
{
    double totlike=0.0;
    //cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            totlike += log(pi(i)*PMAT(i,j))*SEQ2DSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

/*Inference: Calculate the likelihood for divergence t, joint selected nucleotide distribution*/
double likelihoodforTSEQ2DSFS(double t)
{
    return likeSEQwithtwoDSFS(GLOBSEQ2DSFS, t,  GLOBpar);
}

/*Inference: Estimation of divergence time t based on joint selected nucleotide distribution*/
double estimateTSEQ2DSFS(double SEQ2DSFS[4][4], double *t, double parameters[])
{
    /*OK this is stupid*/
    for(int i=0;i<4;i++){
        for (int j=0; j<4; j++){
            GLOBSEQ2DSFS[i][j]=SEQ2DSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ2DSFS, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

//Inference: One EM step for joint selected nucleotide distribution
void EMStepforNuc2DSFS(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites, size_t &psum){
    double p[4][4], psumsr;
//    size_t psum = 0;
    for (int i = 0;i<4;i++){
        for (int j = 0;j<4;j++){
            ESEQSFS2[i][j] = 0;
        }
    }
    for (size_t s = 0; s < numsites; s++){
        for (int r1 = 0;r1<P0[s].size();r1++){
            for (int r2 = 0;r2<P1[s].size();r2++){
                psumsr = 0;
                for (int i = 0;i<4;i++){
                    for (int j=0;j<4;j++){
                        p[i][j] = 0;
                    }
                }
                for (int i = 0;i<4;i++){
                    for (int j=0;j<4;j++){
                        p[i][j] = SEQ2DSFS[i][j]*P0[s][r1].vec[i]*P1[s][r2].vec[j];
                        psumsr +=p[i][j];
                    }
                }
                for (int i = 0;i<4;i++){
                    for (int j=0;j<4;j++){
                        ESEQSFS2[i][j] += p[i][j]/psumsr;
                    }
                }
                psum = psum+1;
            }
        }
    }
//    if (psum > 0){
//        for (int i = 0;i<4;i++){
//            for (int j=0;j<4;j++){
//                ESEQSFS2[i][j] = ESEQSFS2[i][j]/(double)psum;
//            }
//        }
//    }
}

void *bthread(void *aptr){
    EMjobforSEQ2DSFS* ptr = (EMjobforSEQ2DSFS*) aptr;
    size_t index = ptr->index;
    size_t psum = ptr->psum;
    vector<vector<double4> > &p0 = ptr->p0;
    vector<vector<double4> > &p1 = ptr->p1;
    size_t start = ptr -> start;
    size_t len = ptr->len;
    EMStepforNuc2DSFS(ptr->SEQ2DSFS, ptr->segESEQSFS2, p0, p1, len, ptr->psum);
//    ptr->sum = sum;
//    fprintf(stderr,"[%s] thread %lu: Position from %lu to %lu: Finished\n",__FUNCTION__,index,start,start+len-1);
    return NULL;
}

void EMStepforNuc2DSFS_threading_initial(size_t numsites, int nthreads,vector<EMjobforSEQ2DSFS> &jobvec, vector<vector<double4> > &P0, vector<vector<double4> > &P1){
     size_t n = (floor)((double)numsites/(double)nthreads);
     EMjobforSEQ2DSFS tmp;
     tmp.index = 0;
     tmp.psum = 0;
     tmp.start = 0;
     tmp.len = n;
     for (int j=0;j<4;j++){
         for (int k=0;k<4;k++){
             tmp.SEQ2DSFS[j][k] = 0.0625;
         }
     }
    for (size_t i = tmp.start; i < tmp.start+tmp.len; i++){
        tmp.p0.push_back(P0[i]);
        tmp.p1.push_back(P1[i]);
    }
    jobvec.push_back(tmp);
    tmp.p0.clear();
    tmp.p1.clear();
     for (size_t k=1;k<nthreads-1;k++){
         tmp.index = k;
         //tmp.gls =GLDATA;
         //tmp.p_start = testvec.back().p_start+n;
         tmp.start = jobvec.back().start+n;
         tmp.len = n;
         for (size_t i = tmp.start; i < tmp.start+tmp.len; i++){
             tmp.p0.push_back(P0[i]);
             tmp.p1.push_back(P1[i]);
         }
         jobvec.push_back(tmp);
         tmp.p0.clear();
         tmp.p1.clear();
//         cout << k << "\n";
//         if (jobvec[k].p0[0].size()>0 && jobvec[k].p1[0].size()>0){
//         cout<<jobvec[k].p0[0][0].vec[0]<<" "<<jobvec[k].p0[0][0].vec[1]<<" "<<jobvec[k].p0[0][0].vec[2]<<" "<<jobvec[k].p0[0][0].vec[3]<<"\n";
//         cout<<jobvec[k].p1[0][0].vec[0]<<" "<<jobvec[k].p1[0][0].vec[1]<<" "<<jobvec[k].p1[0][0].vec[2]<<" "<<jobvec[k].p1[0][0].vec[3]<<"\n";
//         }
     }
     tmp.index = nthreads-1;
     //tmp.p_start = testvec.back().p_start+n;
     //tmp.gls =GLDATA;
     tmp.start = jobvec.back().start+n;
     tmp.len = numsites-(nthreads-1)*n;
     for (size_t i = tmp.start; i < tmp.start+tmp.len; i++){
        tmp.p0.push_back(P0[i]);
        tmp.p1.push_back(P1[i]);
     }
     jobvec.push_back(tmp);
     tmp.p0.clear();
     tmp.p1.clear();
     P0.clear();
     P1.clear();
     //cout<<testvec.size();
     fprintf(stderr,"[%s] beginning\n",__FUNCTION__);
 }
 
 void EMStepforNuc2DSFS_threading(pthread_t* mythd, int nthreads,vector<EMjobforSEQ2DSFS> &jobvec,double SEQ2DSFS[4][4], size_t numsites){
     size_t i,j,k;
     for(i=0;i<nthreads;i++){
         int rc = pthread_create(&mythd[i],NULL,bthread,&jobvec[i]);
         if(rc)
             fprintf(stderr,"problem starting threadid: %lu\n",i);
     }
     for(i=0;i<nthreads;i++){
         if(pthread_join(mythd[i], NULL)) {
             fprintf(stderr, "Error joining thread\n");
         }
     }
     for (i=0; i<4; i++){
         for (j=0; j<4; j++){
             SEQ2DSFS[i][j]=0.0;
         }
     }
     size_t psum = 0.0;
     for (int i=0; i<nthreads; i++){
         psum += jobvec[i].psum;
         jobvec[i].psum = 0;
     }
     for (j=0;j<4;j++){
         for (k=0;k<4;k++){
             for(i=0;i<nthreads;i++){
                 SEQ2DSFS[j][k] += jobvec[i].segESEQSFS2[j][k];
             }
             SEQ2DSFS[j][k] = SEQ2DSFS[j][k]/(double)psum;
             for(i=0;i<nthreads;i++){
                 jobvec[i].SEQ2DSFS[j][k] = SEQ2DSFS[j][k];
             }
         }
     }
 }


// Inference: Accelerate EM step for joint selected nucleotide distribution
int EMAccelforNuc2DSFS(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites){
    double EM1[4][4], EM2[4][4];
    double stepMin =1;
    double stepMax0 = 1;
    double mstep=4;
    static double stepMax=stepMax0;
    size_t psum1 = 0;
    EMStepforNuc2DSFS(SEQ2DSFS, EM1, P0, P1, numsites, psum1);
    if (psum1 > 0){
        for (int i = 0;i<4;i++){
            for (int j=0;j<4;j++){
                EM1[i][j] = EM1[i][j]/(double)psum1;
            }
        }
    }
    
    double differr1[4][4];
    differr2D(&EM1[0][0], &SEQ2DSFS[0][0], &differr1[0][0], 4);
    double err1 = err2D(&differr1[0][0],4);
    if(sqrt(err1)<tole){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                ESEQSFS2[i][j] = EM1[i][j];
            }
        }
        return 1;
    }
    
    size_t psum2 = 0;
    EMStepforNuc2DSFS(EM1, EM2, P0, P1, numsites, psum2);
    if (psum2 > 0){
        for (int i = 0;i<4;i++){
            for (int j=0;j<4;j++){
                EM2[i][j] = EM2[i][j]/(double)psum2;
            }
        }
    }
    double differr2[4][4];
    differr2D(&EM2[0][0], &EM1[0][0], &differr2[0][0], 4);
    double err2 = err2D(&differr2[0][0],4);
    if(sqrt(err2)<tole){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                ESEQSFS2[i][j] = EM2[i][j];
            }
        }
        return 1;
    }
    
    double differr3[4][4];
    differr2D(&differr2[0][0], &differr1[0][0], &differr3[0][0], 4);
    double err3 = err2D(&differr3[0][0], 4);
    double alpha = sqrt(err1/err3);
    alpha = max(stepMin,min(stepMax,alpha));
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            ESEQSFS2[i][j] = SEQ2DSFS[i][j]+2*alpha*differr1[i][j]+pow(alpha,2)*differr3[i][j];
        }
    }

    int outofparspace =0;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            if(ESEQSFS2[i][j]<0){
                outofparspace++;
                // break;
            }
        }
    }
    if(outofparspace){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                ESEQSFS2[i][j] = EM2[i][j];
            }
        }
    }
    
    if (fabs(alpha - 1) > 0.01){
        //        niter++;
        double ESEQSFS2_tmp[4][4];
        size_t psum_tmp = 0;
        EMStepforNuc2DSFS(ESEQSFS2, ESEQSFS2_tmp, P0, P1, numsites, psum_tmp);
        if (psum_tmp > 0){
            for (int i = 0;i<4;i++){
                for (int j=0;j<4;j++){
                    ESEQSFS2_tmp[i][j] = ESEQSFS2_tmp[i][j]/(double)psum_tmp;
                }
            }
        }
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                swap(ESEQSFS2[i][j],ESEQSFS2_tmp[i][j]);
            }
        }
    }
    //
    if ((alpha - stepMax) > -0.001) {
        stepMax = mstep*stepMax;
    }
    return 1;
}

//Inference: Main EM algorithm for joint genotype distribution
double estimateNuc2DSFS_EM(double SEQ2DSFS[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites)
{
    double ESEQ2DSFS2[4][4], differr[4][4], d;
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            SEQ2DSFS[i][j]=1.0/16.0;
        }
    }
    //tic();
    do {
        size_t psum = 0;
        EMStepforNuc2DSFS(SEQ2DSFS, ESEQ2DSFS2, P0, P1, numsites,psum);
        if (psum > 0){
            for (int i = 0;i<4;i++){
                for (int j=0;j<4;j++){
                    ESEQ2DSFS2[i][j] = ESEQ2DSFS2[i][j]/(double)psum;
                }
            }
        }
        //EMAccelforNuc2DSFS(SEQ2DSFS, ESEQ2DSFS2, P0, P1, numsites);
        differr2D(&ESEQ2DSFS2[0][0], &SEQ2DSFS[0][0], &differr[0][0], 4);
        d = err2D(&differr[0][0],4);
        for (int i = 0; i<4; i++){
            for (int j=0; j<4; j++){
                SEQ2DSFS[i][j] = ESEQ2DSFS2[i][j];
            }
        }
        //printf("d: %lf\n",d);
    }while (sqrt(d)> tole);
    //toc();
    return d;
}

double estimateNuc2DSFS_EM_threading(double SEQ2DSFS[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites,int nthreads)
 {
     double ptemp[4][4], d;
     for (int i=0; i<4; i++){
         for (int j=0; j<4; j++){
             SEQ2DSFS[i][j]=0.0625;
         }
     }
     
     vector<EMjobforSEQ2DSFS> jobvec;
     EMStepforNuc2DSFS_threading_initial(numsites, nthreads, jobvec, P0, P1);
//     EMStepfor2DSFS_threading_initial(numsites,rowL,nthreads,jobvec,GLDATA,SDATA);
     pthread_t *mythd = new pthread_t[nthreads];
//     cout << jobvec.size()<<"\n";
     
     //tic();
     do {
         EMStepforNuc2DSFS_threading(mythd,nthreads,jobvec,ptemp,numsites);
//         for (int i=0; i<4; i++){
//             for (int j=0; j<4; j++){
//                 cout << ptemp[i][j]<<"\t";
//             }
//             cout<<"\n";
//         }
         d=0.0;
         for (int i=0; i<4; i++){
             for (int j=0; j<4; j++){
                 d += pow(ptemp[i][j]-SEQ2DSFS[i][j],2);
                 SEQ2DSFS[i][j]=ptemp[i][j];
             }
         }
         //printf("d: %lf\n",d);
     }while (sqrt(d)> tole);
     
     delete [] mythd;
     //toc();
     return d;
 }

/* Simulation + Inference: Simulation and estimation of divergence time t based on joint selected nucleotide distribution*/
double testsimSEQ2DSFS(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double par[9], int isthreading)
{
    double t, MLV, parameters[8];
    
    
    //SetSeed(6);
    SetSeed(rand()%30000+1);
    
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
//    for (int i=0; i<4; i++)
//    pi(i)=0.25;
//    for (int i=0; i<5; i++)
//    parameters[i]=1.0;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);

//    pi(0)=0.2184;
//    pi(1)=0.2606;
//    pi(2)=0.3265;
//    pi(3)=0.1946;
//    parameters[0]=2.0431;
//    parameters[1]=0.0821;
//    parameters[2]=0.0000;
//    parameters[3]=0.0670;
//    parameters[4]=0.0000;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
    
    
    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    SEQDATA = (double **) malloc(numsites * sizeof(double *));
    //    for (int i = 0; i < numsites; i++)
    //    SEQDATA[i] = (double *) malloc(8 * sizeof(double));
    vector<vector<double4> >P0;
    vector<vector<double4> >P1;
    simulateGLsTwoSpeciesSEQ_v2(RD, numsites, errorrate,  tdiv,  t1,  t2, P0, P1, pijtGTR, parameters);
    
    //Only use for checking.
    size_t effect_numsites = CheckSites(P0, P1, numsites);
    
    double SEQ2DSFS[4][4];
    
    
    
    if (isthreading==1){
        estimateNuc2DSFS_EM_threading(SEQ2DSFS, P0, P1, numsites, 25);
    }else{
        estimateNuc2DSFS_EM(SEQ2DSFS, P0, P1, numsites);
    }
    
    estimateTSEQ2DSFS(SEQ2DSFS, &t, parameters);
    cout<<"Estimated t = "<<t<<"\n";
    
    for (size_t s = 0; s < numsites; s++){
        P0[s].clear();
        P1[s].clear();
    }
    P0.clear();
    P1.clear();
    return t;
}

void GetEMMatrix(double tdiv0, double t10, double t20, size_t numsites, double** Psum, double** GLDATA, int *SDATA){
    for (int i=0;i<40;i++){
        for (int j=0;j<40;j++){
            Psum[i][j]=0;
        }
    }
    gettransitionprobmatGTR(tdiv0);
    Matrix4d PMAT0=PMAT;
    gettransitionprobmatGTR(t10);
    Matrix4d PMAT1=PMAT;
    gettransitionprobmatGTR(t20);
    Matrix4d PMAT2=PMAT;
    for (size_t s=0;s<numsites;s++){
        if (SDATA[s]>0){
            double sum_p=0;
            double Psum0[40][40]={0};
            for (int h1=0;h1<4;h1++){
                double p1 = pi(h1);
                for (int h2=0;h2<4;h2++){
                    double p2 = PMAT0(h1,h2);
                    int g1[2], g2[2];
                    for (int g1_i=0;g1_i<10;g1_i++){
                        findgenotypes_from_index(g1_i, g1);
                        double p3=GLDATA[s][g1_i];
                        if(g1[0]==g1[1]){
                            p3 = p3*pow(PMAT1(h1,g1[0]),2);
                        }else{
                            p3 = p3*2*PMAT1(h1,g1[0])*PMAT1(h1,g1[1]);
                        }
                        for (int g2_j=0;g2_j<10;g2_j++){
                            findgenotypes_from_index(g2_j, g2);
                            //double p = p*GLDATA[s][10+g2_j]*PMAT2(h2,g2[0])*PMAT1(h2,g2[1]);
                            double p=p1*p2*p3*GLDATA[s][10+g2_j];
                            if(g2[0]==g2[1]){
                                p = p*pow(PMAT2(h2,g2[0]),2);
                            }else{
                                p = p*2*PMAT2(h2,g2[0])*PMAT2(h2,g2[1]);
                            }
                            Psum0[h1*10+g1_i][h2*10+g2_j] = Psum0[h1*10+g1_i][h2*10+g2_j]+p;
                            sum_p = sum_p+p;
                        }
                    }
                }
            }
            for (int h1=0;h1<4;h1++){
                for (int h2=0;h2<4;h2++){
                    for (int g1_i=0;g1_i<10;g1_i++){
                        for (int g2_j=0;g2_j<10;g2_j++){
                            Psum[h1*10+g1_i][h2*10+g2_j]=Psum[h1*10+g1_i][h2*10+g2_j]+Psum0[h1*10+g1_i][h2*10+g2_j]/sum_p;
                        }
                    }
                }
            }
            //cout<<s<<"sum_p is "<<sum_p<<"\n";
        }
    }
}


double Q1(double t1, int l){
    gettransitionprobmatGTR(t1);
    double sum0=0;
    if (l==1){
        for (int h1=0;h1<4;h1++){
            int g1[2];
            for (int g1_i=0;g1_i<10;g1_i++){
                findgenotypes_from_index(g1_i, g1);
                double p=0;
                for (int h2=0;h2<4;h2++){
                    for (int g2_j=0;g2_j<10;g2_j++){
                        p = p+GLOBPsum[h1*10+g1_i][h2*10+g2_j];
                    }
                }
                if (g1[0]==g1[1]){
                    sum0 = sum0+(2*log(PMAT(h1,g1[0])))*p;
                }else{
                    sum0 = sum0+(log(PMAT(h1,g1[0]))+log(PMAT(h1,g1[1]))+log(2))*p;
                }
            }
        }
    } else if (l==2){
        for (int h1=0;h1<4;h1++){
            int g1[2];
            for (int g1_i=0;g1_i<10;g1_i++){
                findgenotypes_from_index(g1_i, g1);
                double p=0;
                for (int h2=0;h2<4;h2++){
                    for (int g2_j=0;g2_j<10;g2_j++){
                        p = p+GLOBPsum[h2*10+g2_j][h1*10+g1_i];
                    }
                }
                if (g1[0]==g1[1]){
                    sum0 = sum0+(2*log(PMAT(h1,g1[0])))*p;
                }else{
                    sum0 = sum0+(log(PMAT(h1,g1[0]))+log(PMAT(h1,g1[1]))+log(2))*p;
                }
            }
        }
    } else {
        cout<<"Error: l can be either 1 or 2!";
    }
    return -sum0;
}

double Q11(double t1){
    return Q1(t1,globl1);
}
double Q12(double t1){
    return Q1(t1,globl2);
}

double Q2(double tdiv){
    gettransitionprobmatGTR(tdiv);
    double sum0=0;
    for (int h1=0;h1<4;h1++){
        for (int h2=0;h2<4;h2++){
            double p = 0;
            double ltdiv = log(PMAT(h1,h2));
            for (int g1_i=0;g1_i<10;g1_i++){
                for (int g2_j=0;g2_j<10;g2_j++){
                    p = p+GLOBPsum[h1*10+g1_i][h2*10+g2_j];
                }
            }
            sum0=sum0+ltdiv*p;
        }
    }
    return -sum0;
}

void EMStepforT(double* t0, double *t, size_t numsites,double **GLDATA,int *SDATA){
    GetEMMatrix(t0[0], t0[1], t0[2], numsites, GLOBPsum, GLDATA, SDATA);
    brent(0.0000001, 0.1, 10.0, Q2,  0.000001, &t[0]);
    brent(0.0000001, 0.1, 10.0, Q11, 0.000001, &t[1]);
    brent(0.0000001, 0.1, 10.0, Q12, 0.000001, &t[2]);
}

void differr(double *v1, double *v2, double *v3, int dim){
    for (int i=0;i<dim;i++){
        v3[i] = v1[i]-v2[i];
    }
}

double err(double *v, int dim){
    double e = 0;
    for (int i=0;i<dim;i++){
        e=e+pow(v[i],2);
    }
    return e;
}

int EMAccelforT(double* t0, double *t, size_t numsites,double **GLDATA, int *SDATA){
    double *em1,*em2;
    double stepMin =1;
    double stepMax0 = 1;
    double mstep=4;
    static double stepMax=stepMax0;
    em1 = (double *) malloc(3 * sizeof(double));
    EMStepforT(t0,em1,numsites,GLDATA,SDATA);
    double differr1[3];
    differr(em1, t0, differr1, 3);
    double err1 = err(differr1,3);
    if(sqrt(err1)<tole){
        for(int i=0;i<3;i++){
            t[i]=em1[i];
        }
        return 1;
    }
    
    em2 = (double *) malloc(3 * sizeof(double));
    EMStepforT(em1,em2,numsites,GLDATA,SDATA);
    double differr2[3];
    differr(em2, em1, differr2, 3);
    double err2 = err(differr2,3);
    if(sqrt(err2)<tole){
        for(int i=0;i<3;i++){
            t[i]=em2[i];
        }
        return 1;
    }
    
    double differr3[3];
    differr(differr2, differr1, differr3, 3);
    double err3 = err(differr3,3);
    double alpha = sqrt(err1/err3);
    alpha = std::max(stepMin,std::min(stepMax,alpha));
    for(int i=0;i<3;i++){
        t[i] = t0[i]+2*alpha*differr1[i]+alpha*alpha*differr3[i];
    }
    
    int outofparspace =0;
    for(int i=0;i<3;i++){
        if(t[i]<0){
            outofparspace++;
            // break;
        }
    }
    if(outofparspace){
        for(int i=0;i<3;i++)
        t[i] = em2[i];
    }
    
    if (fabs(alpha - 1) > 0.01){
        //        niter++;
        double *t_tmp;
        t_tmp = (double *) malloc(3 * sizeof(double));
        EMStepforT(t,t_tmp,numsites,GLDATA,SDATA);
        for(int i=0;i<3;i++)
        std::swap(t[i],t_tmp[i]);
    }
    //
    if ((alpha - stepMax) > -0.001) {
        stepMax = mstep*stepMax;
    }
    
    //cout<<t[0]<<" "<<t[1]<<" "<<t[2]<<"\n";
    return 1;
}


void EMlikelihoodforT(double **GLDATA, int* SDATA, size_t numsites, double tdiv0, double t10, double t20, double* t){
    double* t0;
    t0 = (double *) malloc(3 * sizeof(double));
    //tic();
    t0[0]=tdiv0;
    t0[1]=t10;
    t0[2]=t20;
    
    for (int i=0;i<3;i++){
        t[i]=t0[i]+1;
        //cout<<t[i]<<"\n";
    }
    double tdiff[3];
    differr(t, t0, tdiff, 3);
    double tdiff_SumSquare = err(tdiff,3);
    while (sqrt(tdiff_SumSquare)>tole){
        for (int i=0;i<3;i++){
            t0[i] = t[i];
        }
        //EMStepforT(t0, t, numsites,GLDATA,SDATA);
        EMAccelforT(t0, t, numsites,GLDATA,SDATA);
        differr(t, t0, tdiff, 3);
        tdiff_SumSquare = err(tdiff,3);
    }
    cout << t[0] <<" "<<t[1]<<" "<<t[2]<<"\n";
    //toc();
}

double testjointEM(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double* tt){
    
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (size_t i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    GLOBPsum = (double**) malloc(40 * sizeof(double ***));
    for (int i = 0; i < 40; i++){
        GLOBPsum[i] = (double *) malloc(40 * sizeof(double));
    }
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
//    pi(0)=0.2184;
//    pi(1)=0.2606;
//    pi(2)=0.3265;
//    pi(3)=0.1946;
//    parameters[0]=2.0431;
//    parameters[1]=0.0821;
//    parameters[2]=0.0000;
//    parameters[3]=0.0670;
//    parameters[4]=0.0000;
//    for (int i=5;i<8; i++)
//    parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    EMlikelihoodforT(GLDATA, SDATA, numsites, 0.1, 0.04, 0.02,&tt[0]);
    t = tt[0]+tt[1]+tt[2];
    cout<<"Estimated t = "<<t<<"\n";
    
    return t;
}

int gls_writer_double(const char* glfname, int dobinary, int nsites, double** gls){
    //this block opens the filehandle called fp
    
    BGZF *fp = NULL;
//    if(dobinary)
        fp = bgzf_open(glfname,"wb");
//    else
//        fp = bgzf_open(glfname,"wb");

    //for the textoutput we use kstring as internal text buffer
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL;
    kstr->l = kstr->m = 0;
    for(int i=0;i<nsites;i++){
        if(dobinary==1)
            bgzf_write(fp,gls[i],sizeof(double)*20);
        else{
            for(int j=0;j<19;j++){
                ksprintf(kstr,"%f\t",gls[i][j]);
            }
            ksprintf(kstr,"%f\n",gls[i][19]);
            if(kstr->l>10000000){
                bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }

        }
    }
    if(dobinary == 0){
        bgzf_write(fp,kstr->s,kstr->l);
        cout<<"Double txt glffile is constructed!\n";
    }else{
        cout<<"Double binary glffile is constructed!\n";
    }
    bgzf_close(fp);
    return 0;
}

int gls_writer_uchar(const char* glfname, int dobinary, int nsites, uchar **gls){
    //this block opens the filehandle called fp
    BGZF *fp = NULL;
//    if(dobinary)
        fp = bgzf_open(glfname,"wb");
//    else
//        fp = bgzf_open(glfname,"wb");
    
    //for the textoutput we use kstring as internal text buffer
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL;
    kstr->l = kstr->m = 0;
    
    for(int i=0;i<nsites;i++){
        if(dobinary==1)
            bgzf_write(fp,gls[i],sizeof(uchar)*20);
        else{
            for(int j=0;j<19;j++){
                ksprintf(kstr,"%d\t",(int)gls[i][j]);
            }
            ksprintf(kstr,"%d\n",(int)gls[i][19]);
            if(kstr->l>10000000){
                bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            
        }
    }
    if(dobinary == 0){
        bgzf_write(fp,kstr->s,kstr->l);
        cout<<"Unsigned char txt glffile is constructed!\n";
    }else{
        cout<<"Unsigned char binary glffile is constructed!\n";
    }
    bgzf_close(fp);
    return 0;
}

