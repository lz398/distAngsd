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


int main(int argc,char**argv){
    if (argc == 1){
        fprintf(stderr,"\t-> ./distAngsd -o -method -model -glf -vcf -simrep -is2Dinfer -p_inv -isthreading -dobinary -isuchar -numsites -RD -errorrate -tdiv -t1 -t2 -par\n");
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
        int is2Dinfer = p->is2Dinfer;
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
        
        if (p->simrep > 0){
            int simrep = p->simrep;
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
            if (is2Dinfer==1){
                str1 = str1 + "\tis2Dinfer" + "\tp_inv";
                str2 = str2 + "\t" + to_string(is2Dinfer) + "\t" + to_string(p_inv);
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
            if (is2Dinfer==1){
//                if ((!strcasecmp(model,"GTR")) && (strcasecmp(method,"ConsensusGT"))){
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
//                }else{
//                    string repstr = "Symmetric models, such as JC, are proved to have multi-solutions in 2D inferences, please try GTR.\nConsensusGT is not suitable for 2D inferences, please try other methods.\n";
//                    cout << repstr;
//                    if(dobinary)
//                        bgzf_write(fp,repstr.c_str(),repstr.size());
//                    else{
//                        ksprintf(kstr,"%s",repstr.c_str());
//                        bgzf_write(fp,kstr->s,kstr->l);
//                        kstr->l = 0;
//                    }
//                }
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
//            cout<<"Model\tMethod\tVcf\tIshtreading\tIsuchar\n";
//            cout<<model<<"\t"<<method<<"\t"<<vcfname<<"\t"<<isthreading<<"\t"<<isuchar<<"\n";
            string strmodel(model);
            string strmethod(method);
            string strvcfname(vcfname);
            string str1 = "Model\tMethod\tVcf\tIshtreading\tIsuchar";
            string str2 = strmodel+"\t"+strmethod+"\t"+strvcfname+"\t"+to_string(isthreading)+"\t"+to_string(isuchar);
            if (is2Dinfer==1){
                str1 = str1 + "\tis2Dinfer";
                str2 = str2 + "\t" + to_string(is2Dinfer);
            }
            string str = str1+"\n"+str2+"\n";
            cout << str;
            if(dobinary)
                bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            
            
            double t = 0;
            double p = 0;
            string pl=string("PL");
            string fr=string("AFngsrelate");
            char *reg=NULL;
            vcftwoDSFS(vcfname, glfname, 2, 0.00, pl, fr,reg, par, isthreading,isuchar,dobinary,is2Dinfer,t,p);
            string vcfstr;
            if (is2Dinfer == 1){
                cout<<"Estimated t = "<<t<<".\t"<<"Estimated p = "<<p<<".\n";
                vcfstr = "Estimated t = " + to_string(t) + ".\t" + "Estimated p = " + to_string(p) + ".\n";
            }else{
                cout<<"Estimated t = "<<t<<"\n";
                vcfstr = "Estimated t = " + to_string(t) + ".\n";
            }
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
