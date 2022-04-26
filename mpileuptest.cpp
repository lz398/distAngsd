//
//  mpileuptest.cpp
//  
//
//  Created by Lei Zhao on 24/11/2021.
//
#include <iostream>
#include "shared.h"
#include "GLtest.h"
#include "GL2Dtest.h"
#include "vcftest.h"
#include "mpileuptest.h"
using namespace std;
#define LENS 20480

vector<string> split(char* c, const char* src){
    vector<string> res;
    char *saveptr;
    //char* test = strtok_r(c,src);
    //cout << c.find(src) <<"\n";
    char* pch=strchr(c,src[0]);
    int pos=(pch == NULL ? -1 : pch - c);
    if (pos==0){
        res.push_back("");
    }
    char* test = strtok_r(c,src, &saveptr);
    string str0(test);
    res.push_back(str0);
    while(((test = strtok_r(NULL,src, &saveptr)))){
        //while(((test = strtok_r(NULL,"\n\t ")))){
        //cout<<test<<"\t";
        string str1(test);
        res.push_back(str1);
    }
    //cout<<"\n";
    return res;
}

void op(vector<string> &src, const char* opc){
    int sz = src.size();
    if (sz>1){
        if (strcmp(opc,"^")==0){
            for (int i=1;i<sz;i++){
                src[i]=src[i].substr(1);
            }
        }else if(strcmp(opc,"+")==0 || strcmp(opc,"-")==0){
            for (int i=1;i<sz;i++){
                int j = 0;
                while (isdigit(src[i].at(j))!=0){
                    j = j+1;
                }
                src[i]=src[i].substr(j+atoi(src[i].substr(0,j).c_str()));
            }
        }
    }
}

char* join(vector<string> src){
    char* test;
    string res = "";
    if (src.size() == 0){
        test = strdup(res.c_str());
    }else{
        vector<string>::iterator it = src.begin();
        res += *it;
        for (it++;it!=src.end();it++){
            res += *it;
        }
        test = strdup(res.c_str());
    }
    return test;
}

char* readproc(char* c, int &idep, char*& qscores){
    char* mpil;
    //cout<<idep<<"\n";
    if (idep>0){
        vector<string> res0 = split(c, "^");
        op(res0,"^");
        mpil=join(res0);
        vector<string> res1 = split(mpil, "$");
        mpil=join(res1);
        vector<string> res2 = split(mpil, "+");
        op(res2,"+");
        mpil=join(res2);
        vector<string> res3 = split(mpil, "-");
        op(res3,"-");
        mpil=join(res3);
        //cout<<c<<"\n";
        //cout<<mpil<<strlen(mpil)<<"\n";
        for (int n=0;n<strlen(mpil);n++){
            if (mpil[n]-'A'==0 || mpil[n]-'a'==0){
                //cout<<"HaHa"<<"\n";
                //qscores[n]=char(qscores[n]-33);
                mpil[n]=char(33);
                //mpil[n]=char(0);
            }else if(mpil[n]-'C'==0 || mpil[n]-'c'==0){
                //qscores[n]=char(qscores[n]-33);
                mpil[n]=char(34);
            }else if(mpil[n]-'G'==0 || mpil[n]-'g'==0){
                //qscores[n]=char(qscores[n]-33);
                mpil[n]=char(35);
            }else if(mpil[n]-'T'==0 || mpil[n]-'t'==0){
                //qscores[n]=char(qscores[n]-33);
                mpil[n]=char(36);
            }else{
                mpil[n]=char(32);
                qscores[n]=char(32);
                idep--;
            }
        }
        //cout<<idep<<"\n";
        //cout<<(int)mpil[0]<<strlen(mpil)<<"\n";
        if(idep>0){
            string tests(1,char(32));
            vector<string> res4=split(mpil, tests.c_str());
            mpil=join(res4);
            vector<string> res5=split(qscores, tests.c_str());
            qscores=join(res5);
        }
    }else{
        mpil = c;
    }
    return mpil;
}

int mpileupreader(gzFile gz, vector<kstring_t> &indiv0, vector<kstring_t> &indiv1, vector<kstring_t> &bq0, vector<kstring_t> &bq1){
    char buf[LENS];
    char *saveptr0;
    while(gzgets(gz,buf,LENS)){
        if(strlen(buf)>20480-5){
            fprintf(stderr,"\t-> Problem parsing inputfile, increase buffer size\n");
            return 1;
        }
        char *chr = strtok_r(buf,"\n\t ", &saveptr0);
        int pos = atoi(strtok_r(NULL,"\n\t ", &saveptr0));
        char ref = strtok_r(NULL,"\n\t ", &saveptr0)[0];
        //fprintf(stderr,"chr: %s pos: %d ref: %c\n",chr,pos,ref);
        char *tmp = NULL;
        int intid = 0;
        const char* mpil0, *mpil1, *qscores0, *qscores1;
        int k=0;
        int dep0, dep1;
        while(((tmp = strtok_r(NULL,"\n\t ", &saveptr0)))){
            int idep = atoi(tmp);
            char *mpil = strtok_r(NULL,"\n\t ", &saveptr0);
            char *qscores = strtok_r(NULL,"\n\t ", &saveptr0);
            if (k==0){
                //cout<<idep<<" "<<mpil<<" "<<qscores<<"\n";
                mpil0=readproc(mpil,idep,qscores);
                //cout<<mpil<<"\n";
                //cout<<"qscores are "<<qscores<<"\n";
                qscores0 = qscores;
                //cout<<"qscores0 are "<<qscores0<<"\n";
                dep0 = idep;
                k=1;
            }else{
                mpil1=readproc(mpil,idep,qscores);
                qscores1 = qscores;
                dep1 = idep;
                k=0;
                if(dep0>0 && dep1>0){
                    kstring_t read0 = {0,0,NULL};
                    kstring_t q0 = {0,0,NULL};
                    kputs(mpil0,&read0);
                    kputs(qscores0,&q0);
                    kstring_t read1 = {0,0,NULL};
                    kstring_t q1 = {0,0,NULL};
                    kputs(mpil1,&read1);
                    kputs(qscores1,&q1);
                    // Some function here
                    //for (int n=0;n<read0.l;n++){
                    //   if (read0.s[n]-'A'==0 || read0.s[n]-'C'==0 || read0.s[n]-'G'==0 || read0.s[n]-'T' ==0){
                    //       continue;
                    //   }else if(read0.s[n]-'a'==0 || read0.s[n]-'c'==0 || read0.s[n]-'g'==0 || read0.s[n]-'t' ==0){
                    //       read0.s[n]=toupper(read0.s[n]);
                    //   }else{
                    //       read0.s[n]=char(32);
                    //       q0.s[n]=char(32);
                    //   }
                    //}
                    //string tests(1,char(32));
                    //vector<string> res=split(read0.s, tests.c_str());
                    //read0.s=join(res);
                    //cout<<"\n";
                    indiv0.push_back(read0);
                    bq0.push_back(q0);
                    indiv1.push_back(read1);
                    bq1.push_back(q1);
                }
            }
        }
    }
    return 0;
}

double bqplmatrix[PHREDMAX*2][PHREDMAX*2];

void bqplmatrixbuilder(){
    for(int i=0;i<PHREDMAX-1;i++){
        pl2ln[i] = -log(10.0)/(10.0)*i;
        for (int j=0;j<i;j++){
            bqplmatrix[i][j] = (1-exp(pl2ln[i]))*(1-exp(pl2ln[j]));
            bqplmatrix[j][i] = bqplmatrix[i][j];
            bqplmatrix[i+PHREDMAX][j] = (exp(pl2ln[i])/3.0)*(1-exp(pl2ln[j]));
            bqplmatrix[j][i+PHREDMAX] = bqplmatrix[j][i+PHREDMAX];
            bqplmatrix[i][j+PHREDMAX] = (1-exp(pl2ln[i]))*(exp(pl2ln[j])/3.0);
            bqplmatrix[j+PHREDMAX][i] = bqplmatrix[i][j+PHREDMAX];
            bqplmatrix[i+PHREDMAX][j+PHREDMAX] = exp(pl2ln[i]+pl2ln[j])/9.0;
            bqplmatrix[j+PHREDMAX][i+PHREDMAX] = bqplmatrix[i+PHREDMAX][j+PHREDMAX];
        }
        bqplmatrix[i][i] = pow(1-exp(pl2ln[i]),2);
        bqplmatrix[i+PHREDMAX][i] = (exp(pl2ln[i])/3.0)*(1-exp(pl2ln[i]));
        bqplmatrix[i][i+PHREDMAX] = bqplmatrix[i][i+PHREDMAX];
        bqplmatrix[i+PHREDMAX][i+PHREDMAX] = exp(2.0*pl2ln[i])/9.0;
    }
    //pl2ln[PHREDMAX-1] = -log(10.0)/(10.0)*PHREDMAX;
    pl2ln[PHREDMAX-1] = -INFINITY;
    for (int j=0;j<PHREDMAX-1;j++){
        bqplmatrix[PHREDMAX-1][j] = 1-exp(pl2ln[j]);
        bqplmatrix[j][PHREDMAX-1] = bqplmatrix[PHREDMAX-1][j];
        bqplmatrix[2*PHREDMAX-1][j] = 0;
        bqplmatrix[j][2*PHREDMAX-1] = bqplmatrix[2*PHREDMAX-1][j];
        bqplmatrix[PHREDMAX-1][j+PHREDMAX] = exp(pl2ln[j])/3.0;
        bqplmatrix[j+PHREDMAX][PHREDMAX-1] = bqplmatrix[PHREDMAX-1][j+PHREDMAX];
        bqplmatrix[2*PHREDMAX-1][j+PHREDMAX] = 0;
        bqplmatrix[j+PHREDMAX][2*PHREDMAX-1] = bqplmatrix[2*PHREDMAX-1][j+PHREDMAX];
    }
    bqplmatrix[PHREDMAX-1][PHREDMAX-1] = 1;
    bqplmatrix[2*PHREDMAX-1][PHREDMAX-1] = 0;
    bqplmatrix[PHREDMAX-1][2*PHREDMAX-1] = bqplmatrix[2*PHREDMAX-1][PHREDMAX-1];
    bqplmatrix[2*PHREDMAX-1][2*PHREDMAX-1] = bqplmatrix[2*PHREDMAX-1][PHREDMAX-1];
}

int isequal(int a, int b){
    if (a-33==b){
        return 0;
    }else{
        return 1;
    }
}

void EMStepforNuc2DSFS_mpileup(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<kstring_t> *indiv0, vector<kstring_t> *indiv1, vector<kstring_t> *bq0, vector<kstring_t> *bq1, size_t start, size_t numsites, size_t &psum){
    double p[4][4], psumsr;
//    size_t psum = 0;
    for (int i = 0;i<4;i++){
        for (int j = 0;j<4;j++){
            ESEQSFS2[i][j] = 0;
        }
    }
    for (size_t s = start; s < start+numsites; s++){
        for (int r1 = 0;r1<strlen(indiv0->at(s).s);r1++){
            for (int r2 = 0;r2<strlen(indiv1->at(s).s);r2++){
                psumsr = 0;
                for (int i = 0;i<4;i++){
                    for (int j=0;j<4;j++){
                        p[i][j] = 0;
                    }
                }
                for (int i = 0;i<4;i++){
                    for (int j=0;j<4;j++){
                        p[i][j] = SEQ2DSFS[i][j]*bqplmatrix[bq0->at(s).s[r1]-33+isequal(indiv0->at(s).s[r1],i)*PHREDMAX][bq1->at(s).s[r2]-33+isequal(indiv1->at(s).s[r2],j)*PHREDMAX];
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
}

//Inference: Main EM algorithm for joint genotype distribution
double estimateNuc2DSFS_EM_mpileup(double SEQ2DSFS[4][4], vector<kstring_t>* indiv0, vector<kstring_t>* indiv1, vector<kstring_t>* bq0, vector<kstring_t>* bq1, size_t numsites)
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
        EMStepforNuc2DSFS_mpileup(SEQ2DSFS, ESEQ2DSFS2, indiv0, indiv1, bq0, bq1, 0, numsites, psum);
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



// threading is wanted
struct EMjobforSEQ2DSFS_mpileup{
    size_t index;
    size_t psum;
    vector<kstring_t>* indiv0;
    vector<kstring_t>* indiv1;
    vector<kstring_t>* bq0;
    vector<kstring_t>* bq1;
//    vector<> p0;
//    vector<> p1;
    size_t start;
    size_t len;
    double segESEQSFS2[4][4];
    double SEQ2DSFS[4][4];
};

void *dthread(void *aptr){
    EMjobforSEQ2DSFS_mpileup* ptr = (EMjobforSEQ2DSFS_mpileup*) aptr;
//    size_t index = ptr->index;
//    size_t psum = ptr->psum;
//    vector<kstring_t > &indiv0 = ptr->indiv0;
//    vector<kstring_t > &indiv1 = ptr->indiv1;
//    vector<kstring_t > &bp0 = ptr->bp0;
//    vector<kstring_t > &bp1 = ptr->bp1;
//      size_t start = ptr -> start;
//    size_t len = ptr->len;

    EMStepforNuc2DSFS_mpileup(ptr->SEQ2DSFS, ptr->segESEQSFS2, ptr->indiv0, ptr->indiv1, ptr->bq0, ptr->bq1, ptr->start, ptr->len, ptr->psum);
//    ptr->sum = sum;
//    fprintf(stderr,"[%s] thread %lu: Position from %lu to %lu: Finished\n",__FUNCTION__,index,start,start+len-1);
    return NULL;
}
//
void EMStepforNuc2DSFS_mpileup_threading_initial(size_t numsites, int nthreads, vector<EMjobforSEQ2DSFS_mpileup> &jobvec, vector<kstring_t>* indiv0, vector<kstring_t>* indiv1, vector<kstring_t>* bq0, vector<kstring_t>* bq1){
     size_t n = (floor)((double)numsites/(double)nthreads);
     EMjobforSEQ2DSFS_mpileup tmp;
     tmp.index = 0;
     tmp.psum = 0;
     tmp.start = 0;
     tmp.len = n;
     for (int j=0;j<4;j++){
         for (int k=0;k<4;k++){
             tmp.SEQ2DSFS[j][k] = 0.0625;
         }
     }
    tmp.indiv0 = indiv0;
    tmp.indiv1 = indiv1;
    tmp.bq0 = bq0;
    tmp.bq1 = bq1;
//    for (size_t i = tmp.start; i < tmp.start+tmp.len; i++){
//        tmp.p0.push_back(P0[i]);
//        tmp.p1.push_back(P1[i]);
//    }
    jobvec.push_back(tmp);
//    tmp.p0.clear();
//    tmp.p1.clear();
    
     for (size_t k=1;k<nthreads-1;k++){
         tmp.index = k;
         //tmp.gls =GLDATA;
         //tmp.p_start = testvec.back().p_start+n;
         tmp.start = jobvec.back().start+n;
         tmp.len = n;
         jobvec.push_back(tmp);
//         cout << k << "\n";
//         if (jobvec[k].p0[0].size()>0 && jobvec[k].p1[0].size()>0){
//         cout<<jobvec[k].p0[0][0].vec[0]<<" "<<jobvec[k].p0[0][0].vec[1]<<" "<<jobvec[k].p0[0][0].vec[2]<<" "<<jobvec[k].p0[0][0].vec[3]<<"\n";
//         cout<<jobvec[k].p1[0][0].vec[0]<<" "<<jobvec[k].p1[0][0].vec[1]<<" "<<jobvec[k].p1[0][0].vec[2]<<" "<<jobvec[k].p1[0][0].vec[3]<<"\n";
//         }
     }
     tmp.index = nthreads-1;
     tmp.start = jobvec.back().start+n;
     tmp.len = numsites-(nthreads-1)*n;
     jobvec.push_back(tmp);
     //cout<<testvec.size();
     fprintf(stderr,"[%s] beginning\n",__FUNCTION__);
 }

 void EMStepforNuc2DSFS_mpileup_threading(pthread_t* mythd, int nthreads, vector<EMjobforSEQ2DSFS_mpileup> &jobvec, double SEQ2DSFS[4][4], size_t numsites){
     size_t i,j,k;
     for(i=0;i<nthreads;i++){
         int rc = pthread_create(&mythd[i],NULL,dthread,&jobvec[i]);
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

double estimateNuc2DSFS_mpileup_EM_threading(double SEQ2DSFS[4][4], vector<kstring_t>* indiv0, vector<kstring_t>* indiv1, vector<kstring_t>* bq0, vector<kstring_t>* bq1, size_t numsites,int nthreads)
 {
     double ptemp[4][4], d;
     for (int i=0; i<4; i++){
         for (int j=0; j<4; j++){
             SEQ2DSFS[i][j]=0.0625;
         }
     }
     
     vector<EMjobforSEQ2DSFS_mpileup> jobvec;
     EMStepforNuc2DSFS_mpileup_threading_initial(numsites, nthreads, jobvec, indiv0, indiv1, bq0, bq1);
//     EMStepforNuc2DSFS_threading_initial(numsites, nthreads, jobvec, P0, P1);
//     EMStepfor2DSFS_threading_initial(numsites,rowL,nthreads,jobvec,GLDATA,SDATA);
     pthread_t *mythd = new pthread_t[nthreads];
//     cout << jobvec.size()<<"\n";
     
     //tic();
     do {
         EMStepforNuc2DSFS_mpileup_threading(mythd, nthreads, jobvec, ptemp, numsites);
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

int mpileuptwoDSFS(gzFile gz, double par[9], int isthreading, int isuchar, int dobinary, int is2Dinfer, double &t, double &p){
    //size_t nsites = 0;
    vector<kstring_t> indiv0, indiv1, bq0, bq1;
    
    double parameters[8], SEQ2DSFS[4][4];
    t = 0.0;
    p = -1.0;
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    
    int iswrong = mpileupreader(gz, indiv0, indiv1, bq0, bq1);
    if (iswrong !=0){
        return 0;
    }else{
        bqplmatrixbuilder();
        size_t numsites = indiv0.size();
        if (isthreading == 1){
            estimateNuc2DSFS_mpileup_EM_threading(SEQ2DSFS, &indiv0, &indiv1, &bq0, &bq1, numsites, 25);
        }else{
            estimateNuc2DSFS_EM_mpileup(SEQ2DSFS, &indiv0, &indiv1, &bq0, &bq1, numsites);
        }
        if (is2Dinfer == 1){
            double x[2];
            estimateTSEQ2DSFSWithInvSite(SEQ2DSFS, x, parameters);
            cout << "The 2D inferred divergence time is " << x[0] << ",\n";
            cout << "The 2D inferred fraction of invariable sites is " << x[1] << ".\n";
            t = x[0];
            p = x[1];
        }else{
            estimateTSEQ2DSFS(SEQ2DSFS, &t, parameters);
        }
        for (int i=0;i<indiv0.size();i++){
            free(indiv0[i].s);
            free(indiv1[i].s);
            free(bq0[i].s);
            free(bq1[i].s);
        }
        return 0;
    }
}
