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
#include<htslib/bgzf.h>
#include<htslib/kstring.h>
#include "ExistingMethods.h"
#include "io.h"
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

//int main(){
//    
////    for (int i = 0; i<1000; i++){
////        double r = Poisson(0.5);
////        if (r>=0){
////            cout<<r<<"\n";}
////    }
//    //srand (time(NULL));
//    size_t numsites = 10000000;
//    double tdiv =1;
//    double t1 = 0.4;
//    double t2 = 0.25;
//    double errorrate = 0.002;
//    double RD = 1;
//    double p_inv = 0.90;
//    int rep = 1;
//    for (int r = 0; r<rep; r++){
//    cout << "Replicates: " << r<<":\n";
//    SeedSetup();
//    double t, p;
//    testtwoDSFSWithInvSite(RD,numsites, p_inv, tdiv, t1, t2, errorrate, t, p);
//    cout << "The inferred divergence time is " << t << ",\n";
//    cout << "The inferred fraction of invariable sites is " << 100*p << "%.\n";
//
//    double tSEQ, pSEQ;
//    testsimSEQ2DSFSWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ, pSEQ);
//    cout << "The inferred divergence time is " << tSEQ << ",\n";
//    cout << "The inferred fraction of invariable sites is " << 100*pSEQ << "%.\n";
//    
//    
//    double tSEQ1, pSEQ1;
//    testsimSEQDATA_randomWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ1, pSEQ1);
//    cout << "The inferred divergence time is " << tSEQ1 << ",\n";
//    cout << "The inferred fraction of invariable sites is " << 100*pSEQ1 << "%.\n";
//    
//    double tSEQ2, pSEQ2;
//    testsimSEQDATA_consensusWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, tSEQ2, pSEQ2);
//    cout << "The inferred divergence time is " << tSEQ2 << ",\n";
//    cout << "The inferred fraction of invariable sites is " << 100*pSEQ2 << "%.\n";
//    
//    double t_1, p_1;
//    testGL_consensusWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, t_1, p_1);
//    cout << "The inferred divergence time is " << t_1 << ",\n";
//    cout << "The inferred fraction of invariable sites is " << 100*p_1 << "%.\n";
//}
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
//    return 0;
//}

int main(int argc,char**argv){
    if (argc == 1){
        fprintf(stderr,"\t-> ./distAngsd -o -method -model -glf -vcf -simrep -is2Dsim -p_inv -isthreading -dobinary -numsites -RD -errorrate -tdiv -t1 -t2 -par\n");
        fprintf(stderr,"\t-> Default method is geno, default model is JC.\n");
        fprintf(stderr,"\t-> Default number of sites for simulation (numsites) is 1000000, default mean read depth (RD) is 1.0, default calling error rate (errorrate) is 0.002,\n");
        fprintf(stderr,"\t-> default divergent time (tdiv) is 1.0, default t1 is 0.4, and default t2 is 0.25.\n");
        return 0;
    }
    pars *p = get_pars(--argc,++argv);
    if (p!=NULL){
        const char* outname = p->outname;
        const char* method = p->method;
        const char* model = p->model;
        const char* glfname = p->glfname;
        const char* vcfname = p->vcfname;
        int isthreading = p->isthreading;
        int dobinary = p->dobinary;
        double par[9];
        
        BGZF *fp = NULL;
        fp = bgzf_open(outname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        for (int i = 0; i<9; i++){
            par[i]=p->par[i];
            cout << par[i] << " ";
        }
        cout<<"\n";
        if(dobinary)
            bgzf_write(fp,par,sizeof(double)*9);
        else{
            for(int i = 0; i<8; i++){
                ksprintf(kstr,"%f\t",par[i]);
            }
            ksprintf(kstr,"%f\n",par[8]);
            bgzf_write(fp,kstr->s,kstr->l);
            kstr->l = 0;
        }
        
        if (p->simrep!=NULL){
            int simrep = p->simrep;
            int is2Dsim = p->is2Dsim;
            double RD = p->RD;
            int numsites = p->numsites;
            double errorrate = p->errorrate;
            double tdiv = p->tdiv;
            double t1 = p->t1;
            double t2 = p->t2;
            double p_inv = p->p_inv;
            double t = 0;
            double p = 0;
            string str1 = "Model\tMethod\tReplication\tError\ttdiv\tt1\tt2";
            string strmodel(model);
            string strmethod(method);
            string str2 = strmodel+"\t"+strmethod+"\t"+to_string(simrep)+"\t"+to_string(errorrate)+"\t"+to_string(tdiv)+"\t"+to_string(tdiv)+"\t"+to_string(t1)+"\t"+to_string(t2);
            if (is2Dsim==1){
                str1 = str1 + "\tp_inv";
                str2 = str2 + "\t" + to_string(p_inv);
            }
            string str = str1+"\n"+str2+"\n";
            cout<<str;
            if(dobinary)
                bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            if (is2Dsim==1){
                if ((!strcasecmp(model,"GTR")) && (strcasecmp(method,"ConsensusGT"))){
                    for (int r=0; r<simrep; r++){
                        SeedSetup();
                        //string repstr = "Random seed is " + to_string(seed) + "\n";
                        if (!strcasecmp(method,"geno")){
                            testtwoDSFSWithInvSite(RD,numsites, p_inv, tdiv, t1, t2, errorrate, t, p, par, glfname, isthreading, dobinary, r);
                            //t=testtwoDSFS(RD, numsites, tdiv, t1, t2, errorrate, par, glfname, isthreading, dobinary, r);
                        }else if(!strcasecmp(method,"nuc")){
                            testsimSEQ2DSFSWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, t, p, par, isthreading);
                            //t=testsimSEQ2DSFS(RD, numsites, tdiv, t1, t2, errorrate, par);
                        }else if(!strcasecmp(method,"RandomSEQ")){
                            testsimSEQDATA_randomWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, t, p, par);
                            //t=testsimSEQDATA_random(RD, numsites, tdiv, t1, t2, errorrate,par);
                        }else if(!strcasecmp(method,"ConsensusSEQ")){
                            testsimSEQDATA_consensusWithInvSite(RD, numsites, p_inv, tdiv, t1, t2, errorrate, t, p, par);
                            //t=testsimSEQDATA_consensus(RD, numsites, tdiv, t1, t2, errorrate,par);
                        }
                        //                    else if(!strcasecmp(method,"ConsensusGT")){
                        //                        //t=testtwoDSFS_consensusGT(RD, numsites, tdiv, t1, t2, errorrate,par);
                        //                    }
                        //                repstr = repstr + "Estimated t = " + to_string(t) + "\n";
                        string repstr = "Estimated t = " + to_string(t) + ".\t" + "Estimated p = " + to_string(p) + ".\n";
                        if(dobinary)
                            bgzf_write(fp,repstr.c_str(),repstr.size());
                        else{
                            ksprintf(kstr,"%s",repstr.c_str());
                            if(kstr->l>10000000){
                                bgzf_write(fp,kstr->s,kstr->l);
                                kstr->l = 0;
                            }
                        }
                    }
                    if(dobinary == 0){
                        bgzf_write(fp,kstr->s,kstr->l);
                    }
                }else{
                    string repstr = "Symmetric models, such as JC, are proved to have multi-solutions in 2D inferences, please try GTR.\nConsensusGT is not suitable for 2D inferences, please try other methods.\n";
                    cout << repstr;
                    if(dobinary)
                        bgzf_write(fp,repstr.c_str(),repstr.size());
                    else{
                        ksprintf(kstr,"%s",repstr.c_str());
                        bgzf_write(fp,kstr->s,kstr->l);
                        kstr->l = 0;
                    }
                }
            }else{
                for (int r=0; r<simrep; r++){
                    SeedSetup();
                    //string repstr = "Random seed is " + to_string(seed) + "\n";
                    if (!strcasecmp(method,"geno")){
                        t=testtwoDSFS(RD, numsites, tdiv, t1, t2, errorrate, par, glfname, isthreading, dobinary, r);
                    }else if(!strcasecmp(method,"nuc")){
                        t=testsimSEQ2DSFS(RD, numsites, tdiv, t1, t2, errorrate, par, isthreading);
                    }else if(!strcasecmp(method,"RandomSEQ")){
                        t=testsimSEQDATA_random(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }else if(!strcasecmp(method,"ConsensusSEQ")){
                        t=testsimSEQDATA_consensus(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }else if(!strcasecmp(method,"ConsensusGT")){
                        t=testtwoDSFS_consensusGT(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }
                    //                repstr = repstr + "Estimated t = " + to_string(t) + "\n";
                    string repstr = "Estimated t = " + to_string(t) + ".\n";
                    if(dobinary)
                        bgzf_write(fp,repstr.c_str(),repstr.size());
                    else{
                        ksprintf(kstr,"%s",repstr.c_str());
                        if(kstr->l>10000000){
                            bgzf_write(fp,kstr->s,kstr->l);
                            kstr->l = 0;
                        }
                    }
                }
                if(dobinary == 0){
                    bgzf_write(fp,kstr->s,kstr->l);
                }
            }
        }else{
            char* vcfname=p->vcfname;
            int isuchar=p->isuchar;
            cout<<"Model\tMethod\tVcf\tIshtreading\tIsuchar\n";
            cout<<model<<"\t"<<method<<"\t"<<vcfname<<"\t"<<isthreading<<"\t"<<isuchar<<"\n";
            string strmodel(model);
            string strmethod(method);
            string strvcfname(vcfname);
            string str1 = "Model\tMethod\tVcf\tIshtreading\tIsuchar\n";
            string str2 = strmodel+"\t"+strmethod+"\t"+strvcfname+"\t"+to_string(isthreading)+"\t"+to_string(isuchar)+"\n";
            
            string str = str1+str2;
            double t = 0;
            string pl=string("PL");
            string fr=string("AFngsrelate");
            char *reg=NULL;
            t = vcftwoDSFS(vcfname, glfname, 2, 0.00, pl, fr,reg, par, isthreading,isuchar,dobinary);
            cout<<"Estimated t = "<<t<<"\n";
            string vcfstr = "Estimated t = " + to_string(t) + "\n";
            if(dobinary)
                bgzf_write(fp,vcfstr.c_str(),vcfstr.size());
            else{
                ksprintf(kstr,"%s",vcfstr.c_str());
                bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
        }
        bgzf_close(fp);
    }
    return 0;
}
