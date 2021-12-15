#include <iostream>
#include "shared.h"
#include "GLtest.h"
#include "GL2Dtest.h"
#include "bfgs.h"
#include "brent.h"

///* General simulation*/
void simnucleotidesinv(int Genotype[2],double simmat[4][4])
{
    double u;
    int j;
    
    u = uniform();
    j=0;
    while (u>simmat[j][3]) // pi[j]
    {
        j++;
        if (j>3) {printf("Numerical error 1 in simulation algorithm"); exit(-1);}
    }
    Genotype[0]=j;
    Genotype[1]=j;
}

//Simulation: Simulate two individuals' genome and reads in two species.
/*tdiv is the divergence time. t1 and t2 are the average coalescence times within species*/
void simulateGLsTwoSpeciesWithInvSite(double RD, size_t numsites, double p_inv, double errorate, double tdiv, double t1, double t2, double **GLDATA,  double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, *vSITE, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    std::cout << tdiv << " " << t1 << " " << t2 << "\n";
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    vSITE = (int *) malloc(numsites*(sizeof(int)));
    for (size_t i=0; i<numsites; i++)
    ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    int ** TrueGenome = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (size_t i=0; i<numsites; i++)
    TrueGenome[i]=(int *) malloc(4*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    double u;
//    int nuccount[4];
//    for (int i = 0; i < 4; i++){
//        nuccount[i] = 0;
//    }
    for (size_t i=0; i<numsites; i++){
        u = uniform();
        if (u<=p_inv){
            simnucleotidesinv(ancDATA[i], SIMMAT);
            //std::cout << ancDATA[i][0] << "\n";
            //nuccount[ancDATA[i][0]] = nuccount[ancDATA[i][0]] + 1;
            vSITE[i] = 0;
        }else{
            simnucleotides(ancDATA[i], SIMMAT);
            vSITE[i] = 1;
        }
        // printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    std::cout<<"Ancient Data was derived!\n";
    //std::cout<<nuccount[0]<<" "<<nuccount[1]<<" "<<nuccount[2]<<" "<<nuccount[3]<<"\n";
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    //std::cout << "Individual 0:\n";
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    //makeGLswitchmatrix(errorate, switchmatrix, simswitchmatrix);
    for (size_t i=0; i<numsites; i++){
        if (vSITE[i] == 1){
            simpoly(ancDATA[i][0], SIMMAT, genotypes);
        }else{
            genotypes[0] = ancDATA[i][0];
            genotypes[1] = ancDATA[i][0];
        }
        //std::cout<<"Ancient nucleotide is "<<ancDATA[i][0]<<"\n";
        //std::cout<<genotypes[0]<<" "<<genotypes[1]<<"\n";
        TrueGenome[i][0] = genotypes[0];
        TrueGenome[i][1] = genotypes[1];
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        //    simGLs(GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix);
        //std::cout<<"Error rate is "<<errorate<<"\n";
        simGLs_reads(GLDATA, genotypes, i, 0, errorate, RD);
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    //std::cout << "Individual 1:\n";
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (size_t i=0; i<numsites; i++){
        if (vSITE[i]==1){
            simpoly(ancDATA[i][1], SIMMAT, genotypes);
        }else{
            genotypes[0] = ancDATA[i][1];
            genotypes[1] = ancDATA[i][1];
        }
        //std::cout<<"Ancient nucleotide is "<<ancDATA[i][1]<<"\n";
        //std::cout<<genotypes[0]<<" "<<genotypes[1]<<"\n";
        TrueGenome[i][2] = genotypes[0];
        TrueGenome[i][3] = genotypes[1];
        //     simGLs(GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix);
        //std::cout<<"Error rate is "<<errorate<<"\n";
        simGLs_reads(GLDATA, genotypes, i, 1, errorate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    /*then we free memory*/
    for (size_t i=0; i<numsites; i++){
    free(ancDATA[i]);
    free(TrueGenome[i]);
    }
    free(ancDATA);
    free(vSITE);
    free(TrueGenome);
}

////Inference: total likelihood based on the inferred joint genotype distribution.
double likeGLwithtwoDSFSWithInvSite(double twoDSFS[10][10], const double* x, double par[8])
{
    int g1[2], g2[2];
    double like, totlike=0.0;
    double t = x[0];
    double p_inv = x[1];
//    std::cout<<"t is "<<t<<"\n";
//    std::cout<<"p_inv is "<<p_inv<<"\n";
    //    std::cout << "t is " << x[0]<<"\n";
    //    std::cout << "p_inv is "<<x[1]<<"\n";
    //std::cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            like =0.0;
            //if (g1[0]!=g2[0] || g1[1]!=g2[1] || g1[0]!=g2[1] || g1[1]!=g2[0] || g1[0]!=g1[1] || g2[0]!=g2[1]){
                for (int k=0; k<2; k++){
                    for (int v=0; v<2; v++){
                        if (g1[k]!=g2[v]){
                            like = like + log(pi(g1[k])*(1-p_inv)*PMAT(g1[k],g2[v]));
                        }else{
                            like = like + log(pi(g1[k])*(p_inv+(1-p_inv)*PMAT(g1[k],g2[v])));
                        }
                    }
                }
            //}else{
            //    like = 4*log(pi(g1[0])*(p_inv+(1-p_inv)*PMAT(g1[0],g2[0])));
            //}
            totlike += like*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

// Partial likelihood
double parlikeGLwithtwoDSFSWithInvSite(double twoDSFS[10][10], double t, double par[8])
{
    int g1[2], g2[2];
    double like, parlike=0.0;
//    std::cout<<"t is "<<t<<"\n";
//    std::cout<<"p_inv is "<<p_inv<<"\n";
    //    std::cout << "t is " << x[0]<<"\n";
    //    std::cout << "p_inv is "<<x[1]<<"\n";
    //std::cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            like =0.0;
            //if (g1[0]!=g2[0] || g1[1]!=g2[1] || g1[0]!=g2[1] || g1[1]!=g2[0] || g1[0]!=g1[1] || g2[0]!=g2[1]){
                for (int k=0; k<2; k++){
                    for (int v=0; v<2; v++){
                        if (g1[k]!=g2[v]){
                            like = like + log(pi(g1[k])*PMAT(g1[k],g2[v]));
                        }
                    }
                }
            //}else{
            //    like = 4*log(pi(g1[0])*(p_inv+(1-p_inv)*PMAT(g1[0],g2[0])));
            //}
            parlike += like*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -parlike;
}



//
void gettransitionprobmatGTR_grad(double t)
{
    VectorXd EXPOSdt(4);
    int i, j, k;
    
    for (k=0; k<4; k++)
    EXPOSdt(k) = RRVAL(k)*exp(t*RRVAL(k));
    for (i=0; i<4; i++){
        for (j=0; j<4; j++){
            PMATdt(i,j) = 0.0;
            for (k=0; k<4; k++)
            PMATdt(i,j) =  PMATdt(i,j) + RRVEC(k,j)*LRVEC(i,k)*EXPOSdt(k);
        }
    }
    //newt=0;
    /*printf("PMAT:\n");
     for (i=0; i<4;i++){
     for (j=0; j<4; j++)
     printf("%lf ",PMAT[i][j]);
     printf("\n");
     }*/
}
//

void likeGLwithtwoDSFSWithInvSite_grad(double twoDSFS[10][10], const double* x, double *y, double par[8]){
    int g1[2], g2[2];
    double t = x[0];
    double p_inv = x[1];
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    gettransitionprobmatGTR_grad(t);
    
    y[0] =0.0;
    y[1] =0.0;
    double a, b;
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            a = 0;
            b = 0;
            //if (g1[0]!=g2[0] || g1[1]!=g2[1] || g1[0]!=g2[1] || g1[1]!=g2[0] || g1[0]!=g1[1] || g2[0]!=g2[1]){
            for (int k=0; k<2; k++){
                for (int v=0; v<2; v++){
                    if (g1[k]!=g2[v]){
                        a = a + PMATdt(g1[k],g2[v])/PMAT(g1[k],g2[v]);
                        b = b - 1/(1-p_inv);
                    }else{
                        a = a + (1-p_inv)*PMATdt(g1[k],g2[v])/(p_inv+(1-p_inv)*PMAT(g1[k],g2[v]));
                        b = b + (1-PMAT(g1[k],g2[v]))/(p_inv+(1-p_inv)*PMAT(g1[k],g2[v]));
                    }
                }
            }
            y[0] += a*twoDSFS[i][j];
            y[1] += b*twoDSFS[i][j];
        }
    }
    y[0] = -y[0];
    y[1] = -y[1];
}

////Inference: Calculate the likelihood for divergence t, joint genotype distribution
int ncalls;
double likelihoodforTWithInvSite(const double* x,const void *)
{
    ncalls++;
    return likeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS, x,  GLOBpar);
}

int ncalls_grad;
void likelihoodforTWithInvSite_grad(const double* x, double* y)
{
    ncalls_grad++;
    likeGLwithtwoDSFSWithInvSite_grad(GLOBtwoDSFS, x, y,  GLOBpar);
}


double GLOB_p_inv;
double likelihoodforTWithInvSitePinvKnown(double t)
{
    double x[2];
    x[0] = t;
    x[1] = GLOB_p_inv;
    return likeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS, x,  GLOBpar);
}

double parlikelihoodforTWithInvSite(double t){
    return parlikeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS, t,  GLOBpar);
}
//
////Inference: Calculate the likelihood for divergence t, joint genotype distribution + tree structure
//double likelihoodforT_m(double t)
//{
//    return likeGLwithtwoDSFS_m(GLOBtwoDSFS, t,  GLOBpar, globt1, globt2);
//}
//
////Inference: Estimation of divergence time t based on joint genotype distribution
double estimateTWithInvSite(double twoDSFS[10][10], double x[2], double parameters[])
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
    
    //define lims
    double lbd[2] = {0.0000001,0.0000001};
    double ubd[2] = {9.9999999,0.9999999};
    int nbd[2] = {2,2};
    //no grad
    double invec[2] = {1.5,0.5};
    
//    ncalls=0;
//    double MLV = findmax_bfgs(2,invec,NULL,likelihoodforTWithInvSite,NULL,lbd,ubd,nbd,-1);
//    x[0] = invec[0];
//    x[1] = invec[1];
//    std::cout << "nfunctioncalls: " << ncalls << "\n";
//    std::cout << "Inferred MLV is " << MLV << "\n";
    
//    double z[2];
    ncalls=0;
    ncalls_grad=0;
    double MLV = findmax_bfgs(2,invec,NULL,likelihoodforTWithInvSite,likelihoodforTWithInvSite_grad,lbd,ubd,nbd,-1);
    x[0] = invec[0];
    x[1] = invec[1];
    std::cout << "nfunctioncalls: " << ncalls << "\n";
    std::cout << "nfunctioncalls_grad: " << ncalls_grad << "\n";
    //std::cout << "Inferred z is " << x[0] << " " << x[1] << "\n";
    std::cout << "Inferred MLV of Geno2DSFS is " << MLV << "\n";
    //fprintf(stderr,"nograd val: %f=(%f,%f) nfunctioncalls: %d\n",nograd,invec[0],invec[1],ncalls);
    //double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

double estimateTWithInvSitePinvKnown(double twoDSFS[10][10], double *t, double parameters[],double p_inv)
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
    
    GLOB_p_inv = p_inv;
  
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTWithInvSitePinvKnown, 0.000001, t);
//    std::cout << "nfunctioncalls: " << ncalls << "\n";
//    std::cout << "nfunctioncalls_grad: " << ncalls_grad << "\n";
    //std::cout << "Inferred z is " << x[0] << " " << x[1] << "\n";
    std::cout << "Inferred MLV of Geno2DSFS with p_inv value "<< p_inv <<" is " << -MLV << "\n";
    //fprintf(stderr,"nograd val: %f=(%f,%f) nfunctioncalls: %d\n",nograd,invec[0],invec[1],ncalls);
    //double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

double estimateTWithInvSiteParlike(double twoDSFS[10][10], double *t, double parameters[]){
    /*OK this is stupid*/
    for(int i=0;i<10;i++){
        for (int j=0; j<10; j++){
            GLOBtwoDSFS[i][j]=twoDSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    double MLV = brent(0.0000001, 0.1, 10.0, parlikelihoodforTWithInvSite, 0.000001, t);
    return MLV;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution*/
void testtwoDSFSWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9], const char* glfname, int isthreading, int dobinary, int r)
{
    double **GLDATA, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (size_t i = 0; i < numsites; i++)
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
 
    simulateGLsTwoSpeciesWithInvSite(RD, numsites, p_inv, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    if (glfname!=NULL){
        string str = glfname+to_string(r);
        gls_writer_double(str.c_str(), dobinary, numsites, GLDATA);
    }
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    //
    
    //Filter the effective numsites
    size_t eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
//    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
//        for (int i = 0; i<10; i++){
//            for (int j = 0; j<10; j++){
//                std::cout << twoDSFS[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//    double p_inv_hat = 0;
    
//    estimate2DSFS_EM_threading(twoDSFS, GLDATA, SDATA, numsites, eff_numsites, 25, 20);
    if (isthreading==1){
        estimate2DSFS_EM_threading(twoDSFS, GLDATA, SDATA, numsites, eff_numsites, 25, 2*10);
    }else{
        estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    }

    
    //std::cout << "The inferred p_inv_hat is " << p_inv_hat <<"\n";
    //Estimate T
    double x[2];
    estimateTWithInvSite(twoDSFS, x, parameters);
    std::cout << "The 2D inferred divergence time is " << x[0] << ",\n";
    std::cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";

    t = x[0];
    p = x[1];


    for (size_t i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);

    free(SDATA);
}


void testGL_consensusWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p)
{
    double **GLDATA, parameters[8], twoDSFS[10][10];
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
    
    //simulate data
    std::cout << "p_inv is " << p_inv << "\n";
    simulateGLsTwoSpeciesWithInvSite(RD, numsites, p_inv, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    size_t eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    int m0, m1, r;
    double u;
    for (int i=0; i<10; i++){
        for(int j=0; j<10; j++){
            twoDSFS[i][j]=0.0;
        }
    }
    for (size_t i=0; i<numsites; i++){
        if (SDATA[i]==1){
            vector<int> n0;
            vector<int> n1;
            double GLtemp0[10], GLtemp1[10];
            for (int j = 0; j<10; j++){
                GLtemp0[j] = GLDATA[i][j];
                GLtemp1[j] = GLDATA[i][j+10];
            }
            findmaxvec(GLtemp0,10,n0);
            findmaxvec(GLtemp1,10,n1);
            
            if (n0.size()==1){
                m0 = n0[0];
            }else{
                u = uniform();
                r = floor(u*n0.size());
                m0 = n0[r];
            }
            if (n1.size()==1){
                m1 = n1[0];
            }else{
                u = uniform();
                r = floor(u*n1.size());
                m1 = n1[r];
            }
        
            twoDSFS[m0][m1] += 1.0;
        }
    }

    
    double p_inv_hat = 0.0;
    for (int i = 0; i<10; i++){
        for (int j = 0; j<10; j++){
            twoDSFS[i][j] = twoDSFS[i][j]/(double)eff_numsites;
            std::cout << twoDSFS[i][j] << "\t";
        }
        if (i==0 || i==4 || i==7 || i==9){
            p_inv_hat += twoDSFS[i][i];
        }
        std::cout << "\n";
    }
    
    //Estimate T
    double x[2];
    estimateTWithInvSite(twoDSFS, x, parameters);
    std::cout << "The 2D inferred divergence time is " << x[0] << ",\n";
    std::cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";
    estimateTWithInvSiteParlike(twoDSFS, &t, parameters);
    std::cout<< "The inferred t via partial likelihood (off-diagnal elements) is "  << t <<".\n";
    estimateTWithInvSitePinvKnown(twoDSFS,&t,parameters,p_inv);
    std::cout<< "The inferred t given p_inv = " << p_inv << " is " << t <<"\n";
    estimateTWithInvSitePinvKnown(twoDSFS,&t,parameters,p_inv_hat);
    std::cout<< "The inferred t given p_inv_hat = " << p_inv_hat << " is " << t <<".\n";
    estimateTWithInvSitePinvKnown(twoDSFS,&t,parameters,x[1]);
    std::cout<< "The inferred t given 2D inferred p_inv = " << x[1] << " is " << t <<".\n";
    const double y[2] = {tdiv,p_inv};
    std::cout<<"True value likelihood is "<<-likeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS, y,  GLOBpar)<<".\n";
//    const double truetest[2] = {tdiv,p_inv};
//    std::cout <<"First " << -likeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS,truetest,GLOBpar)<<"\n";
//    const double truetest1[2] = {x[0],x[1]};
//    std::cout <<"Second " <<  -likeGLwithtwoDSFSWithInvSite(GLOBtwoDSFS,truetest1,GLOBpar)<<"\n";
//    std::cout<<"Estimated t = "<<x[0]<<"\n";
//    std::cout<<"Estimated p_inv = "<<x[1]<<"\n";
    
    t = x[0];
    p = x[1];


    for (size_t i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);

    free(SDATA);
}

void findmaxvec(double vec[], int n, vector<int> &indices){
    double current_max = -1.0;
    for (int i = 0; i < n; i++){
        if (vec[i] > current_max)
        {
            current_max = vec[i];
            indices.clear();
        }
        if (vec[i] == current_max)
        {
            indices.push_back(i);
        }
    }
}
///*Simulation+Inference: Simulation and estimation of divergence time t based on a chosen read per site*/
void testsimSEQDATA_randomWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9])
{
    double MLV, parameters[8], SEQ2DSFS[4][4];


    //SetSeed(6);
    SetSeed(rand()%30000+1);
    
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }

    vector<vector<double4> >P0;
    vector<vector<double4> >P1;
    simulateGLsTwoSpeciesSEQ2DSFSWithInvSite(RD, numsites, p_inv, errorrate, tdiv, t1, t2, P0, P1, pijtGTR, parameters);

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            SEQ2DSFS[i][j] = 0;
        }
    }
    
    size_t eff_numsites=0;
    int m0, m1;
    for (size_t i=0; i<numsites; i++){
        if (P0[i].size()>0 && P1[i].size()>0){
            double u0 = uniform();
            double u1 = uniform();
            int r0 = floor(u0*P0[i].size());
            int r1 = floor(u1*P1[i].size());
            vector<int> n0;
            vector<int> n1;
            findmaxvec(P0[i][r0].vec,4,n0);
            findmaxvec(P1[i][r1].vec,4,n1);
            double u;
            int r, m0, m1;
            if (n0.size()==1){
                m0 = n0[0];
            }else{
                u = uniform();
                r = floor(u*n0.size());
                m0 = n0[r];
            }
            if (n1.size()==1){
                m1 = n1[0];
            }else{
                u = uniform();
                r = floor(u*n1.size());
                m1 = n1[r];
            }
            SEQ2DSFS[m0][m1] += 1.0;
            eff_numsites += 1;
        }
    }
    
//    double p_inv_hat=0;
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 4; j++){
//            SEQ2DSFS[i][j] = SEQ2DSFS[i][j]/(double)eff_numsites;
//            std::cout << SEQ2DSFS[i][j] << "\t";
//        }
//        p_inv_hat += SEQ2DSFS[i][i];
//        std::cout << "\n";
//    }
//    std::cout << "The inferred p_inv_hat is " << p_inv_hat <<"\n";

    double x[2];
    estimateTSEQ2DSFSWithInvSite(SEQ2DSFS, x, parameters);
    std::cout << "The 2D inferred divergence time is " << x[0] << ",\n";
    std::cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";
//    estimateTSEQ2DSFSWithInvSiteParlike(SEQ2DSFS, &t, parameters);
//    std::cout<< "The inferred t via partial likelihood (off-diagnal elements) is "<< t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv);
//    std::cout<< "The inferred t given p_inv = " << p_inv << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv_hat);
//    std::cout<< "The inferred t given p_inv_hat = " << p_inv_hat << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,x[1]);
//    std::cout<< "The inferred t given 2D inferred p_inv = " << x[1] << " is " << t <<".\n";
    t = x[0];
    p = x[1];
//    const double y[2] = {tdiv,p_inv};
//    std::cout<<"True value likelihood is "<<-likeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, y,  GLOBpar)<<".\n";
    
    P0.clear();
    P1.clear();
//    for (int i=0; i<8; i++)
//    GLOBpar[i]=parameters[i];
//    globnumsites=numsites;
//    globerror=errorrate;
//    MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ_v1, 0.000001, &t);
//    std::cout<<"Estimated t = "<<t<<"\n";
//
//    for (int i = 0; i < numsites; i++)
//    free(SEQDATA[i]);
//    free(SEQDATA);
//    free(SEQ_SDATA);
}


void testsimSEQDATA_consensusWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9])
{
    double MLV, parameters[8], SEQ2DSFS[4][4];


    //SetSeed(6);
    SetSeed(rand()%30000+1);
    
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }

    vector<vector<double4> >P0;
    vector<vector<double4> >P1;
    simulateGLsTwoSpeciesSEQ2DSFSWithInvSite(RD, numsites, p_inv, errorrate, tdiv, t1, t2, P0, P1, pijtGTR, parameters);

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            SEQ2DSFS[i][j] = 0;
        }
    }

    size_t eff_numsites=0;
    double v0[4], v1[4], u;
    int m0, m1, mm0, mm1, r;
    for (size_t i=0; i<numsites; i++){
        if (P0[i].size()>0 && P1[i].size()>0){
            for (int k = 0; k<4; k++){
                v0[k] = 0.0;
                v1[k] = 0.0;
            }
            for (int k = 0; k<P0[i].size(); k++){
                vector<int> n0;
                findmaxvec(P0[i][k].vec,4,n0);
                if (n0.size()==1){
                    m0 = n0[0];
                }else{
                    u = uniform();
                    r = floor(u*n0.size());
                    m0 = n0[r];
                }
                v0[m0] -= 10/log(10)*log(1-P0[i][k].vec[m0]);
            }
            for (int k = 0; k<P1[i].size(); k++){
                vector<int> n1;
                findmaxvec(P1[i][k].vec,4,n1);
                if (n1.size()==1){
                    m1 = n1[0];
                }else{
                    u = uniform();
                    r = floor(u*n1.size());
                    m1 = n1[r];
                }
                v1[m1] -= 10/log(10)*log(1-P1[i][k].vec[m1]);
            }
            vector<int> nn0;
            vector<int> nn1;
            findmaxvec(v0,4,nn0);
            findmaxvec(v1,4,nn1);
            if (nn0.size()==1){
                mm0 = nn0[0];
            }else{
                u = uniform();
                r = floor(u*nn0.size());
                mm0 = nn0[r];
            }
//
//            std::cout<<mm0<<"\n";
            if (nn1.size()==1){
                mm1 = nn1[0];
            }else{
                u = uniform();
                r = floor(u*nn1.size());
                mm1 = nn1[r];
            }
            SEQ2DSFS[mm0][mm1] += 1.0;
            eff_numsites += 1;
        }
    }

//    double p_inv_hat=0;
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 4; j++){
//            SEQ2DSFS[i][j] = SEQ2DSFS[i][j]/(double)eff_numsites;
//            std::cout << SEQ2DSFS[i][j] << "\t";
//        }
//        p_inv_hat += SEQ2DSFS[i][i];
//        std::cout << "\n";
//    }
//    std::cout << "The inferred p_inv_hat is " << p_inv_hat <<"\n";


    double x[2];
    estimateTSEQ2DSFSWithInvSite(SEQ2DSFS, x, parameters);
    std::cout << "The 2D inferred divergence time is " << x[0] << ",\n";
    std::cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";
//    estimateTSEQ2DSFSWithInvSiteParlike(SEQ2DSFS, &t, parameters);
//    std::cout<< "The inferred t via partial likelihood (off-diagnal elements) is "<< t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv);
//    std::cout<< "The inferred t given p_inv = " << p_inv << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv_hat);
//    std::cout<< "The inferred t given p_inv_hat = " << p_inv_hat << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,x[1]);
//    std::cout<< "The inferred t given 2D inferred p_inv = " << x[1] << " is " << t <<".\n";
    t = x[0];
    p = x[1];
//    const double y[2] = {tdiv,p_inv};
//    std::cout<<"True value likelihood is "<<-likeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, y,  GLOBpar)<<".\n";

    P0.clear();
    P1.clear();
//    for (int i=0; i<8; i++)
//    GLOBpar[i]=parameters[i];
//    globnumsites=numsites;
//    globerror=errorrate;
//    MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ_v1, 0.000001, &t);
//    std::cout<<"Estimated t = "<<t<<"\n";
//
//    for (int i = 0; i < numsites; i++)
//    free(SEQDATA[i]);
//    free(SEQDATA);
//    free(SEQ_SDATA);
}

/*Simulation: Simulate two individuals' genome and reads in two species.
 tdiv is the divergence time. t1 and t2 are the average coalescence times within species, based on the joint selected nucleotide distribution*/
void simulateGLsTwoSpeciesSEQ2DSFSWithInvSite(double RD, size_t numsites, double p_inv, double errorrate, double tdiv, double t1, double t2, vector<vector<double4> >&P0, vector<vector<double4> >&P1, double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10], u;
    int **ancDATA, *vSITE, genotypes[2];

    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}

    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    vSITE = (int *) malloc(numsites*(sizeof(int)));
    for (size_t i=0; i<numsites; i++)
    ancDATA[i]=(int *) malloc(2*(sizeof(int)));

    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    
    // A little may need to be changed here!
    for (size_t i=0; i<numsites; i++){
        u = uniform();
        if (u<=p_inv){
            simnucleotidesinv(ancDATA[i], SIMMAT);
            //std::cout << ancDATA[i][0] << "\n";
            //nuccount[ancDATA[i][0]] = nuccount[ancDATA[i][0]] + 1;
            vSITE[i] = 0;
        }else{
            simnucleotides(ancDATA[i], SIMMAT);
            vSITE[i] = 1;
        }
        //printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    std::cout<<"Ancient Data was derived!\n";

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
        if (vSITE[i] == 1){
            simpoly(ancDATA[i][0], SIMMAT, genotypes);
        }else{
            genotypes[0] = ancDATA[i][0];
            genotypes[1] = ancDATA[i][0];
        }
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
        if (vSITE[i] == 1){
            simpoly(ancDATA[i][1], SIMMAT, genotypes);
        }else{
            genotypes[0] = ancDATA[i][1];
            genotypes[1] = ancDATA[i][1];
        }
        //simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_v2(P1, genotypes, i, 1, errorrate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    /*then we free memory*/
    for (size_t i=0; i<numsites; i++)
    free(ancDATA[i]);
    free(ancDATA);
    free(vSITE);
}
//
///*Check the effective number of sites*/
//int CheckSites(vector<vector<double4> >&P0, vector<vector<double4> >&P1, int numsites){
//    int eff_numsites = 0;
//    for (int s=0;s<numsites;s++){
//        if (P0[s].size()>0 && P1[s].size()>0){
//            eff_numsites += 1;
//        }
//    }
//    std::cout << "Effective number of sites is "<<eff_numsites<<"\n";
//    return eff_numsites;
//}
//
///*Inference: total likelihood based on the inferred joint selected nucleotide distribution.*/
double likeSEQwithtwoDSFSWithInvSite(double SEQ2DSFS[4][4], const double* x, double par[8])
{
    double totlike=0.0;
    double t = x[0];
    double p_inv = x[1];
    //std::cout<<t<<"\n";

    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/

    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            if (i!=j){
                totlike += log(pi(i)*(1-p_inv)*PMAT(i,j))*SEQ2DSFS[i][j];
            }else{
                totlike += log(pi(i)*(p_inv+(1-p_inv)*PMAT(i,j)))*SEQ2DSFS[i][j];
            }
        }
    }
    
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

double parlikeSEQwithtwoDSFSWithInvSite(double SEQ2DSFS[4][4], const double t, double par[8])
{
    double parlike=0.0;
    //std::cout<<t<<"\n";

    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/

    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            if (i!=j){
                parlike += log(pi(i)*PMAT(i,j))*SEQ2DSFS[i][j];
            }
        }
    }
    
    //    printf("%lf: %lf\n",t,totlike);
    return -parlike;
}

void likeSEQwithtwoDSFSWithInvSite_grad(double SEQ2DSFS[4][4], const double* x, double *y, double par[8]){
    double totlike=0.0;
    double t = x[0];
    double p_inv = x[1];
    //std::cout<<t<<"\n";

    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    gettransitionprobmatGTR_grad(t);
    
    y[0] = 0.0;
    y[1] = 0.0;
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            if (i!=j){
                y[0] += PMATdt(i,j)/PMAT(i,j)*SEQ2DSFS[i][j];
                y[1] -= 1/(1-p_inv)*SEQ2DSFS[i][j];
            }else{
                y[0] += (1-p_inv)*PMATdt(i,j)/(p_inv+(1-p_inv)*PMAT(i,j))*SEQ2DSFS[i][j];
                y[1] += (1-PMAT(i,j))/(p_inv+(1-p_inv)*PMAT(i,j))*SEQ2DSFS[i][j];
            }
        }
    }
    y[0] = -y[0];
    y[1] = -y[1];
}

///*Inference: Calculate the likelihood for divergence t, joint selected nucleotide distribution*/
int ncallsSEQ2DSFS;
double likelihoodforTSEQ2DSFSWithInvSite(const double* x, const void *)
{
    ncallsSEQ2DSFS++;
    return likeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, x,  GLOBpar);
}

int ncallsSEQ2DSFS_grad;
void likelihoodforTSEQ2DSFSWithInvSite_grad(const double* x, double* y)
{
    ncallsSEQ2DSFS_grad++;
    likeSEQwithtwoDSFSWithInvSite_grad(GLOBSEQ2DSFS, x, y,  GLOBpar);
}

double likelihoodforTSEQ2DSFSWithInvSitePinvKnown(double t)
{
    double x[2];
    x[0] = t;
    x[1] = GLOB_p_inv;
    return likeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, x,  GLOBpar);
}

double parlikelihoodforTSEQ2DSFSWithInvSite(double t)
{
    return parlikeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, t, GLOBpar);
}

/*Inference: Estimation of divergence time t based on joint selected nucleotide distribution*/
double estimateTSEQ2DSFSWithInvSite(double SEQ2DSFS[4][4], double x[2], double parameters[])
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

    //define lims
    double lbd[2] = {0.0000001,0.0000001};
    double ubd[2] = {9.9999999,0.9999999};
    int nbd[2] = {2,2};
    //no grad
    double invec[2] = {1.5,0.5};
    
//    ncallsSEQ2DSFS=0;
//    double MLV = findmax_bfgs(2,invec,NULL,likelihoodforTSEQ2DSFSWithInvSite,NULL,lbd,ubd,nbd,-1);
//    x[0] = invec[0];
//    x[1] = invec[1];
//    std::cout << "nfunctioncalls: " << ncallsSEQ2DSFS << "\n";
//    std::cout << "Inferred MLV is " << MLV << "\n";
//    std::cout << "(t,p_inv)= "<<"("<<x[0]<<","<<x[1]<<")\n";

    ncallsSEQ2DSFS=0;
    ncallsSEQ2DSFS_grad=0;
    double MLV = findmax_bfgs(2,invec,NULL,likelihoodforTSEQ2DSFSWithInvSite,likelihoodforTSEQ2DSFSWithInvSite_grad,lbd,ubd,nbd,-1);
    x[0] = invec[0];
    x[1] = invec[1];
    std::cout << "nfunctioncalls: " << ncallsSEQ2DSFS << "\n";
    std::cout << "nfunctioncalls_grad: " << ncallsSEQ2DSFS_grad << "\n";
    //std::cout << "Inferred z is " << x[0] << " " << x[1] << "\n";
    std::cout << "Inferred MLV of 2D Nuc2DSFS is " << MLV << "\n";
//    std::cout << "(t,p_inv)= "<<"("<<x[0]<<","<<x[1]<<")\n";
    return MLV;
}

double estimateTSEQ2DSFSWithInvSitePinvKnown(double SEQ2DSFS[4][4], double* t, double parameters[],double p_inv)
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
    
    GLOB_p_inv = p_inv;
    
//    //define lims
//    double lbd[2] = {0.0000001,0.0000001};
//    double ubd[2] = {9.9999999,0.9999999};
//    int nbd[2] = {2,2};
//    //no grad
//    double invec[2] = {1.5,0.5};
    
//    ncalls=0;
//    double MLV = findmax_bfgs(2,invec,NULL,likelihoodforTWithInvSite,NULL,lbd,ubd,nbd,-1);
//    x[0] = invec[0];
//    x[1] = invec[1];
//    std::cout << "nfunctioncalls: " << ncalls << "\n";
//    std::cout << "Inferred MLV is " << MLV << "\n";
    
//    double z[2];
//    ncalls=0;
//    ncalls_grad=0;
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ2DSFSWithInvSitePinvKnown, 0.000001, t);
//    std::cout << "nfunctioncalls: " << ncalls << "\n";
//    std::cout << "nfunctioncalls_grad: " << ncalls_grad << "\n";
    //std::cout << "Inferred z is " << x[0] << " " << x[1] << "\n";
    std::cout << "Inferred MLV of Nuc2DSFS with p_inv value "<< p_inv <<" is " << -MLV << "\n";
    //fprintf(stderr,"nograd val: %f=(%f,%f) nfunctioncalls: %d\n",nograd,invec[0],invec[1],ncalls);
    //double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

double estimateTSEQ2DSFSWithInvSiteParlike(double SEQ2DSFS[4][4], double* t, double parameters[])
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
    
    double MLV = brent(0.0000001, 0.1, 10.0, parlikelihoodforTSEQ2DSFSWithInvSite, 0.000001, t);
//    std::cout << "nfunctioncalls: " << ncalls << "\n";
//    std::cout << "nfunctioncalls_grad: " << ncalls_grad << "\n";
    //std::cout << "Inferred z is " << x[0] << " " << x[1] << "\n";
//    std::cout << "Inferred MLV of Nuc2DSFS is " << MLV << "\n";
    //fprintf(stderr,"nograd val: %f=(%f,%f) nfunctioncalls: %d\n",nograd,invec[0],invec[1],ncalls);
    //double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}


//
/* Simulation + Inference: Simulation and estimation of divergence time t based on joint selected nucleotide distribution*/
void testsimSEQ2DSFSWithInvSite(double RD, int numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9], int isthreading)
{
    double MLV, parameters[8], x[2];


    //SetSeed(6);
    SetSeed(rand()%30000+1);

    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }


    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    SEQDATA = (double **) malloc(numsites * sizeof(double *));
    //    for (int i = 0; i < numsites; i++)
    //    SEQDATA[i] = (double *) malloc(8 * sizeof(double));
    vector<vector<double4> >P0;
    vector<vector<double4> >P1;
    simulateGLsTwoSpeciesSEQ2DSFSWithInvSite(RD, numsites, p_inv, errorrate, tdiv, t1, t2, P0, P1, pijtGTR, parameters);

    //Only use for checking.
    size_t effect_numsites = CheckSites(P0, P1, numsites);

    double SEQ2DSFS[4][4];
//    double trueSEQ2DSFS[4][4];
//
//    std::cout << "The true SEQ2DSFS is " << "\n";
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 4; j++){
//            if (i == j){
//                trueSEQ2DSFS[i][j] = pi(i)*p_inv;
//            }else{
//                trueSEQ2DSFS[i][j] = 0;
//            }
//            trueSEQ2DSFS[i][j] += pi(i)*PMAT(i,j)*(1-p_inv);
//            std::cout << trueSEQ2DSFS[i][j] << "\t";
//        }
//        std::cout << "\n";
//    }
//    estimateTSEQ2DSFSWithInvSite(trueSEQ2DSFS, x, parameters);
//    std::cout << "True estimators for t and p_inv is " << x[0]<< " "<<x[1]<<"\n";
//    estimateNuc2DSFS_EM(SEQ2DSFS, P0, P1, numsites);
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 4; j++){
//            std::cout << SEQ2DSFS[i][j] << "\t";
//            SEQ2DSFS[i][j] = 0;
//        }
//        std::cout << "\n";
//    }

    
    if (isthreading==1){
        estimateNuc2DSFS_EM_threading(SEQ2DSFS, P0, P1, numsites, 25);
    }else{
        estimateNuc2DSFS_EM(SEQ2DSFS, P0, P1, numsites);
    }
//    double p_inv_hat=0;
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 4; j++){
//            std::cout << SEQ2DSFS[i][j] << "\t";
//        }
//        p_inv_hat += SEQ2DSFS[i][i];
//        std::cout << "\n";
//    }
//    std::cout << "The inferred p_inv_hat is " << p_inv_hat <<"\n";

    estimateTSEQ2DSFSWithInvSite(SEQ2DSFS, x, parameters);
    std::cout << "The 2D inferred divergence time is " << x[0] << ",\n";
    std::cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";
//    estimateTSEQ2DSFSWithInvSiteParlike(SEQ2DSFS, &t, parameters);
//    std::cout<< "The inferred t via partial likelihood (off-diagnal elements) is "<< t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv);
//    std::cout<< "The inferred t given the true p_inv = " << p_inv << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,p_inv_hat);
//    std::cout<< "The inferred t given p_inv_hat = " << p_inv_hat << " is " << t <<".\n";
//    estimateTSEQ2DSFSWithInvSitePinvKnown(SEQ2DSFS, &t, parameters,x[1]);
//    std::cout<< "The inferred t given 2D inferred p_inv = " << x[1] << " is " << t <<".\n";
//    const double y[2] = {tdiv,p_inv};
//    std::cout<<"True value likelihood is "<<-likeSEQwithtwoDSFSWithInvSite(GLOBSEQ2DSFS, y,  GLOBpar)<<".\n";
    
    t = x[0];
    p = x[1];
    
//    for (int s = 0; s < numsites; s++){
//        P0[s].clear();
//        P1[s].clear();
//    }
//    P0.clear();
//    P1.clear();
}
