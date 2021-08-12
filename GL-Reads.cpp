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
#include "shared.h"
#include "vcftest.h"
#include "GLtest.h"
#include "GL2Dtest.h"
#include "ExistingMethods.h"
using namespace std;
using namespace Eigen;

//Inference
double tole=1e-6;

MatrixXd RRVEC(4,4), LRVEC(4,4), PMAT(4,4), PMATdt(4,4), PMATdt2(4,4);
VectorXd RRVAL(4), pi(4); /*RRVEC and LRVEC are globals used by
                           diagonalizeGTR and gettransitionprobmatGTR.
                           PMAT and pi is the transition probability matrix and stationary distcb., respectively, accessed by diagonalizeGTR and gettransitionprobmatGTR and pijtGTR.
                           */
double GLOBtwoDSFS[10][10], GLOBSEQ2DSFS[4][4], GLOBpar[8], globerror, globt1, globt2;
size_t globnumsites;
double **SEQDATA, **GLOBPsum;;
int *SEQ_SDATA;
int **DATA; /*the data is stored here*/
int newt, newmat; /*these are control variables used to determine if the transition probability matrix needs to be recalculated*/

int main(){
    
//    for (int i = 0; i<1000; i++){
//        double r = Poisson(0.5);
//        if (r>=0){
//            cout<<r<<"\n";}
//    }
    //srand (time(NULL));
    size_t numsites = 10000000;
    double tdiv =1;
    double t1 = 0.4;
    double t2 = 0.25;
    double errorrate = 0.002;
    double RD = 1;
    double p_inv = 0.90;
    int rep = 1;
    for (int r = 0; r<rep; r++){
    cout << "Replicates: " << r<<":\n";
    SeedSetup();
    double t, p;
    testtwoDSFSWithInvSite(RD,numsites, p_inv, tdiv, t1, t2, errorrate, t, p);
    cout << "The inferred divergence time is " << t << ",\n";
    cout << "The inferred fraction of invariable sites is " << 100*p << "%.\n";

    double tSEQ, pSEQ;
    testsimSEQ2DSFSWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ, pSEQ);
    cout << "The inferred divergence time is " << tSEQ << ",\n";
    cout << "The inferred fraction of invariable sites is " << 100*pSEQ << "%.\n";
    
    
    double tSEQ1, pSEQ1;
    testsimSEQDATA_randomWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ1, pSEQ1);
    cout << "The inferred divergence time is " << tSEQ1 << ",\n";
    cout << "The inferred fraction of invariable sites is " << 100*pSEQ1 << "%.\n";
    
    double tSEQ2, pSEQ2;
    testsimSEQDATA_consensusWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ2, pSEQ2);
    cout << "The inferred divergence time is " << tSEQ2 << ",\n";
    cout << "The inferred fraction of invariable sites is " << 100*pSEQ2 << "%.\n";
    
    double t_1, p_1;
    testGL_consensusWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, t_1, p_1);
    cout << "The inferred divergence time is " << t_1 << ",\n";
    cout << "The inferred fraction of invariable sites is " << 100*p_1 << "%.\n";
}
//    double m = 0.0;
//    double v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t = testtwoDSFS(RD,numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean1 = m/(double)rep;
//    double var1 = v/(double)rep - pow(mean1,2);
//    cout<<"GL Mean is "<<mean1<<", "<<"Variance is "<<var1<<"\n";
////
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t = testtwoDSFS_m(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean2 = m/(double)rep;
//    double var2 = v/(double)rep - pow(mean2,2);
//    cout<<"GLm Mean is "<<mean2<<", "<<"Variance is "<<var2<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    double mm[3]={0.0,0.0,0.0};
//    double vv[3]={0.0,0.0,0.0};
//    double tt[3];
//    for (int r = 0;r<rep;r++){
//        t = testjointEM(RD, numsites, tdiv, t1, t2, errorrate, &tt[0]);
//        m = m + t;
//        v = v + pow(t,2);
//        for (int i = 0;i<3;i++){
//            mm[i] = mm[i] + tt[i];
//            vv[i] = vv[i] + pow(tt[i],2);
//        }
//    }
//    double mean3 = m/(double)rep;
//    double var3 = v/(double)rep - pow(mean3,2);
//    double mmean3[3];
//    double vvar3[3];
//    for (int i = 0;i<3;i++){
//        mmean3[i] = mm[i]/(double)rep;
//        vvar3[i] = vv[i]/(double)rep - pow(mmean3[i],2);
//    }
//    cout<<"EM Mean is "<<mean3<<", "<<"Variance is "<<var3<<"\n";
//    cout<<"EM Means are "<<mmean3[0]<<" "<<mmean3[1]<<" "<<mmean3[2]<<"\n";
//    cout<<"Variance Means are "<<vvar3[0]<<" "<<vvar3[1]<<" "<<vvar3[2]<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t = testsimSEQDATA_v1(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean4 = m/(double)rep;
//    double var4 = v/(double)rep - pow(mean4,2);
//    cout<<"SEQ Mean is "<<mean4<<", "<<"Variance is "<<var4<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testsimSEQ2DSFS( RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean5 = m/(double)rep;
//    double var5 = v/(double)rep - pow(mean5,2);
//    cout<<"SEQ2DSFS Mean is "<<mean5<<", "<<"Variance is "<<var5<<"\n";
//
//    cout<<"The followings are the comparisons of preexisting methods!\n";
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testtwoDSFS_Ambiguity( RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean6 = m/(double)rep;
//    double var6 = v/(double)rep - pow(mean6,2);
//    cout<<"Ambiguity inference Mean is "<<mean6<<", "<<"Variance is "<<var6<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testtwoDSFS_Ambiguity_v2(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean7 = m/(double)rep;
//    double var7 = v/(double)rep - pow(mean7,2);
//    cout<<"Ambiguity inference Mean is "<<mean7<<", "<<"Variance is "<<var7<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testsimSEQDATA_random(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean8 = m/(double)rep;
//    double var8 = v/(double)rep - pow(mean8,2);
//    cout<<"SEQ_random inference Mean is "<<mean8<<", "<<"Variance is "<<var8<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testsimSEQDATA_consensus(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean9 = m/(double)rep;
//    double var9 = v/(double)rep - pow(mean9,2);
//    cout<<"SEQ_consensus inference Mean is "<<mean9<<", "<<"Variance is "<<var9<<"\n";
//
//    m = 0.0;
//    v = 0.0;
//    for (int r = 0;r<rep;r++){
//        t=testtwoDSFS_consensusGT(RD, numsites, tdiv, t1, t2, errorrate);
//        m = m + t;
//        v = v + pow(t,2);
//    }
//    double mean10 = m/(double)rep;
//    double var10 = v/(double)rep - pow(mean10,2);
//    cout<<"ConsensusGT inference Mean is "<<mean10<<", "<<"Variance is "<<var10<<"\n";
////    for (int i=0;i<28;i++){
////        double tdiv1 = tdiv + (i-1)*0.05;
////        for (int r=0;r<100;r++){
////            testtwoDSFS_ngsDist(RD, numsites, tdiv1, t1, t2, errorrate);
////            testtwoDSFS_ngsDist_snp(RD, numsites, tdiv1, t1, t2, errorrate);
////        }
////    }
    return 0;
}
