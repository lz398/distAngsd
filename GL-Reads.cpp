#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include<htslib/bgzf.h>
#include<htslib/kstring.h>

#include "shared.h"
#include "vcftest.h"
#include "mpileuptest.h"
#include "GLtest.h"
#include "GL2Dtest.h"
#include "ExistingMethods.h"
#include "io.h"

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


std::string to_string_int(int a){
  char buf[1024];
  snprintf(buf,1024,"%d",a);
  return std::string(buf);
}

std::string to_string_float(double d){
  char buf[1024];
  snprintf(buf,1024,"%f",d);
  return std::string(buf);
}


int main(int argc,char**argv){
    if (argc == 1){
        fprintf(stderr,"\t-> ./distAngsd -o -method -model -inglf -outglf -vcf -mpileup -simrep -is2Dinfer -isex -p_inv -isthreading -inbin -inuchar -outbin -outuchar -numsites -RD -errorrate -tdiv -t1 -t2 -par\n");
        fprintf(stderr,"\t-> Default method is geno, default model is JC, default setting considers both transversions and transitions (isex=0).\n");
        fprintf(stderr,"\t-> Default number of sites for simulation (numsites) is 1000000, default mean read depth (RD) is 1.0, default calling error rate (errorrate) is 0.002,\n");
        fprintf(stderr,"\t-> default divergent time (tdiv) is 1.0, default t1 is 0.4, and default t2 is 0.25.\n");
        return 0;
    }
    pars *ptr = get_pars(--argc,++argv);
    if (ptr!=NULL){
        const char* outname = ptr->outname;
        const char* method = ptr->method;
        const char* model = ptr->model;
        const char* glfname = ptr->glfname;
        char* vcfname = ptr->vcfname;
        char* mpileupname = ptr->mpileupname;
        const char* tabname = ptr->tabname;
        int isthreading = ptr->isthreading;
        int dobinary = ptr->dobinary;
        int is2Dinfer = ptr->is2Dinfer;
        int isex = ptr->isex;
        int simrep = ptr->simrep;
        double par[9];
        
        BGZF *fp = NULL;
        fp = bgzf_open(outname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        for (int i = 0; i<9; i++){
            par[i]=ptr->par[i];
            cout << par[i] << " ";
        }
        cout<<"\n";
        if(dobinary)
            my_bgzf_write(fp,par,sizeof(double)*9);
        else{
            for(int i = 0; i<8; i++){
                ksprintf(kstr,"%f\t",par[i]);
            }
            ksprintf(kstr,"%f\n",par[8]);
            my_bgzf_write(fp,kstr->s,kstr->l);
            kstr->l = 0;
        }
        
        if (simrep > 0){
            double RD = ptr->RD;
            int numsites = ptr->numsites;
            double errorrate = ptr->errorrate;
            //cout<<errorrate<<"Hello!\n";
            double tdiv = ptr->tdiv;
            double t1 = ptr->t1;
            double t2 = ptr->t2;
            double p_inv = ptr->p_inv;
            double t = 0;
            double p = 0;
            char str1[4096] = "Model\tMethod\tReplication\tError\ttdiv\tt1\tt2";
	    //            string str2 = strmodel+"\t"+strmethod+"\t"+to_string(simrep)+"\t"+to_string(errorrate)+"\t"+to_string(tdiv)+"\t"+to_string(t1)+"\t"+to_string(t2);

	    char str2[4096];
	    snprintf(str2,4096,"%s\t%s\t%d\t%f\t%f\t%f\t%f",model,method,simrep,errorrate,tdiv,t1,t2);
            if (is2Dinfer==1){
	      //                str1 = str1 + "\tis2Dinfer" + "\tp_inv";
	      //                str2 = str2 + "\t" + to_string(is2Dinfer) + "\t" + to_string(p_inv);
		snprintf(str1,4096,"%s\tis2Dinfer\tp_inv",str1);
		snprintf(str2,4096,"%s\t%d\t%f",str2,is2Dinfer,p_inv);
            }
	    //            str1 = str1 + "\tThreading" + "\tOut_uchar" + "\tOut_binary";
	    //            str2 = str2 +  "\t" + to_string(isthreading) + "\t" + to_string(0) + "\t" + to_string(dobinary);
	    snprintf(str1,4096,"%s\tThreading\tExclude\tOut_uchar\tOut_binary",str1);
	    snprintf(str2,4096,"%s\t%d\t%d\t0\t%d",str2,isthreading,isex,dobinary);
	    
            string str = std::string(str1)+"\n"+std::string(str2)+"\n";
            cout<<str;
            if(dobinary)
                my_bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            if (is2Dinfer==1){
//                if ((!strcasecmp(model,"GTR")) && (strcasecmp(method,"AmbiguityGT"))){
                    for (int r=0; r<simrep; r++){
                        SeedSetup();
                        //string repstr = "Random seed is " + to_string(seed) + "\n";
                        if (!strcasecmp(method,"geno")){
                            testtwoDSFSWithInvSite(RD,numsites, p_inv, tdiv, t1, t2, errorrate, t, p, par, glfname, isthreading, dobinary,r);
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
                        //                    else if(!strcasecmp(method,"AmbiguityGT")){
                        //                        //t=testtwoDSFS_AmbiguityGT(RD, numsites, tdiv, t1, t2, errorrate,par);
                        //                    }
                        //                repstr = repstr + "Estimated t = " + to_string(t) + "\n";
                        string repstr = "Estimated t = " + to_string_float(t) + ".\t" + "Estimated p = " +to_string_float(p) + ".\n";
                        if(dobinary)
                            my_bgzf_write(fp,repstr.c_str(),repstr.size());
                        else{
                            ksprintf(kstr,"%s",repstr.c_str());
                            if(kstr->l>10000000){
                                my_bgzf_write(fp,kstr->s,kstr->l);
                                kstr->l = 0;
                            }
                        }
                    }
                    if(dobinary == 0){
                        my_bgzf_write(fp,kstr->s,kstr->l);
                    }
//                }else{
//                    string repstr = "Symmetric models, such as JC, are proved to have multi-solutions in 2D inferences, please try GTR.\nAmbiguityGT is not suitable for 2D inferences, please try other methods.\n";
//                    cout << repstr;
//                    if(dobinary)
//                        my_bgzf_write(fp,repstr.c_str(),repstr.size());
//                    else{
//                        ksprintf(kstr,"%s",repstr.c_str());
//                        my_bgzf_write(fp,kstr->s,kstr->l);
//                        kstr->l = 0;
//                    }
//                }
            }else{
                for (int r=0; r<simrep; r++){
                    SeedSetup();
                    //string repstr = "Random seed is " + to_string(seed) + "\n";
                    if (!strcasecmp(method,"geno")){
                        t=testtwoDSFS(RD, numsites, tdiv, t1, t2, errorrate, par, glfname, isthreading, dobinary,isex, r);
                    }else if(!strcasecmp(method,"nuc")){
                        t=testsimSEQ2DSFS(RD, numsites, tdiv, t1, t2, errorrate, par, isthreading, isex);
                    }else if(!strcasecmp(method,"RandomSEQ")){
                        t=testsimSEQDATA_random(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }else if(!strcasecmp(method,"ConsensusSEQ")){
                        t=testsimSEQDATA_consensus(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }else if(!strcasecmp(method,"AmbiguityGT")){
                        t=testtwoDSFS_AmbiguityGT(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }else if(!strcasecmp(method,"NoAmbiguityGT")){
                        t=testtwoDSFS_NoAmbiguityGT(RD, numsites, tdiv, t1, t2, errorrate,par);
                    }
                    //                repstr = repstr + "Estimated t = " + to_string(t) + "\n";
                    string repstr = "Estimated t = " + to_string_float(t) + ".\n";
                    if(dobinary)
                        my_bgzf_write(fp,repstr.c_str(),repstr.size());
                    else{
                        ksprintf(kstr,"%s",repstr.c_str());
                        if(kstr->l>10000000){
                            my_bgzf_write(fp,kstr->s,kstr->l);
                            kstr->l = 0;
                        }
                    }
                }
                if(dobinary == 0){
                    my_bgzf_write(fp,kstr->s,kstr->l);
                }
            }
        }else if(vcfname !=NULL){
            int isuchar=ptr->isuchar;
//            cout<<"Model\tMethod\tVcf\tIshtreading\tIsuchar\n";
//            cout<<model<<"\t"<<method<<"\t"<<vcfname<<"\t"<<isthreading<<"\t"<<isuchar<<"\n";
            string strmodel(model);
            string strmethod(method);
            string strvcfname(vcfname);
            string str1 = "Model\tMethod\tVcf";
            string str2 = strmodel+"\t"+strmethod+"\t"+strvcfname;
            if (is2Dinfer==1){
                str1 = str1 + "\tis2Dinfer";
                str2 = str2 + "\t" + to_string_int(is2Dinfer);
            }
            str1 = str1 + "\tThreading" + "\tExclude" + "\tOut_uchar" + "\tOut_binary";
            str2 = str2 +  "\t" + to_string_int(isthreading) + "\t"+ to_string_int(isex) + "\t" + to_string_int(0) + "\t" + to_string_int(dobinary);
            string str = str1+"\n"+str2+"\n";
            cout << str;
            if(dobinary)
                my_bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            
            double t = 0;
            double p = 0;
            string pl=string("PL");
            string fr=string("AFngsrelate");
            char *reg=NULL;
            vcftwoDSFS(vcfname, glfname, 2, 0.00, pl, fr,reg, par, isthreading,isuchar,dobinary,is2Dinfer,isex,t,p);
            string vcfstr;
            if (is2Dinfer == 1){
                cout<<"Estimated t = "<<t<<".\t"<<"Estimated p = "<<p<<".\n";
                vcfstr = "Estimated t = " + to_string_float(t) + ".\t" + "Estimated p = " + to_string_float(p) + ".\n";
            }else{
                cout<<"Estimated t = "<<t<<".\n";
                vcfstr = "Estimated t = " + to_string_float(t) + ".\n";
            }
            if(dobinary)
                my_bgzf_write(fp,vcfstr.c_str(),vcfstr.size());
            else{
                ksprintf(kstr,"%s",vcfstr.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
        }else if(mpileupname != NULL){
            int isuchar=ptr->isuchar;
            string strmodel(model);
            string strmethod(method);
            string strmpileupname(mpileupname);
            string str1 = "Model\tMethod\tMpileup";
            string str2 = strmodel+"\t"+strmethod+"\t"+strmpileupname;
            if (is2Dinfer==1){
                str1 = str1 + "\tis2Dinfer";
                str2 = str2 + "\t" + to_string_int(is2Dinfer);
            }
            str1 = str1 + "\tThreading" + "\tExclude" + "\tOut_uchar" + "\tOut_binary";
            str2 = str2 +  "\t" + to_string_int(isthreading) + "\t" + to_string_int(isex) + "\t" + to_string_int(0) + "\t" + to_string_int(dobinary);
            string str = str1+"\n"+str2+"\n";
            cout << str;
            if(dobinary)
                my_bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            
            gzFile gz = Z_NULL;
            gz = gzopen(mpileupname,"rb");
            double t = 0;
            double p = 0;
            if(gz==Z_NULL){
                fprintf(stderr,"\t-> Problem opening inputfile or inputstream\n");
                return 0;
            }
            mpileuptwoDSFS(gz, par, isthreading, isuchar, dobinary, is2Dinfer, isex, t, p);
            string mpileupstr;
            if (is2Dinfer == 1){
                cout<<"Estimated t = "<<t<<".\t"<<"Estimated p = "<<p<<".\n";
                mpileupstr = "Estimated t = " + to_string_float(t) + ".\t" + "Estimated p = " + to_string_float(p) + ".\n";
            }else{
                cout<<"Estimated t = "<<t<<"\n";
                mpileupstr = "Estimated t = " + to_string_float(t) + ".\n";
            }
            if(dobinary)
                my_bgzf_write(fp,mpileupstr.c_str(),mpileupstr.size());
            else{
                ksprintf(kstr,"%s",mpileupstr.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
        }else if(tabname !=NULL){
            int isuchar=ptr->isuchar;
            int tabbinary=ptr->tabbinary;
            int tabuchar=ptr->tabuchar;
            double t = 0;
            double p = 0;
            
            string strmodel(model);
            string strmethod(method);
            string strtabname(tabname);
            string str1 = "Model\tMethod\tInGLTab\tIn_uchar\tIn_binary";
            string str2 = strmodel+"\t"+strmethod+"\t"+strtabname+"\t"+to_string_int(tabuchar)+"\t"+to_string_int(tabbinary);
            if (is2Dinfer==1){
                str1 = str1 + "\tis2Dinfer";
                str2 = str2 + "\t" + to_string_int(is2Dinfer);
            }
            str1 = str1 + "\tThreading" + "\tExclude" + "\tOut_uchar" + "\tOut_binary";
            str2 = str2 +  "\t" + to_string_int(isthreading) + "\t" + to_string_int(isex)  + "\t" + to_string_int(0) + "\t" + to_string_int(dobinary);
            string str = str1+"\n"+str2+"\n";
            cout << str;
            if(dobinary)
                my_bgzf_write(fp,str.c_str(),str.size());
            else{
                ksprintf(kstr,"%s",str.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
            
            tabletwoDSFS(tabname, par, isthreading, tabuchar, tabbinary, is2Dinfer, isex, t, p);
            string vcfstr;
            if (is2Dinfer == 1){
                cout<<"Estimated t = "<<t<<".\t"<<"Estimated p = "<<p<<".\n";
                vcfstr = "Estimated t = " + to_string_float(t) + ".\t" + "Estimated p = " + to_string_float(p) + ".\n";
            }else{
                cout<<"Estimated t = "<<t<<"\n";
                vcfstr = "Estimated t = " + to_string_float(t) + ".\n";
            }
            if(dobinary)
                my_bgzf_write(fp,vcfstr.c_str(),vcfstr.size());
            else{
                ksprintf(kstr,"%s",vcfstr.c_str());
                my_bgzf_write(fp,kstr->s,kstr->l);
                kstr->l = 0;
            }
        }
        bgzf_close(fp);
    }
    return 0;
}
