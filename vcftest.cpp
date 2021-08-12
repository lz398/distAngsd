#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <cmath>
#include <limits>
#include <string>
#include <pthread.h>
#include <cassert>
#include "shared.h"
#include "GLtest.h"
#include "vcftest.h"
//#include "special.h"
#define diskio_threads 24
using namespace std;

int std_queue = 0;
pthread_mutex_t mymut = PTHREAD_MUTEX_INITIALIZER;
int mycounter = 0; //my semaphore
//typedef struct satan_t{
//    char*fname;
//    int minind;
//    double minfreq;
//    string vcf_format_field;
//    string vcf_allele_field;
//    char *seek;
//    vector<double *> mygl;
//    vector<double> freqs;
//    int nind;
//}satan;

vector<satan> jobs;

int findnuctypeindex(char* c1){
    string s(c1);
    int i = s.compare("A");
    if (i == 0){
        return 0;
    }else if(i == 2){
        return 1;
    }else if(i == 6){
        return 2;
    }else if(i == 19){
        return 3;
    }else if(s.compare("<*>") == 0){
        return 4;
    }
}
//
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
//    }else{printf("error in genotype index conversion table"); exit(-1);}
//
//}
//
//void findgenotypes_from_index(int inde, int genotype[2])
//{
//    switch(inde) {
//        case 0: genotype[0]=0; genotype[1]=0; break;
//        case 1: genotype[0]=0; genotype[1]=1; break;
//        case 2: genotype[0]=0; genotype[1]=2; break;
//        case 3: genotype[0]=0; genotype[1]=3; break;
//        case 4: genotype[0]=1; genotype[1]=1; break;
//        case 5: genotype[0]=1; genotype[1]=2; break;
//        case 6: genotype[0]=1; genotype[1]=3; break;
//        case 7: genotype[0]=2; genotype[1]=2; break;
//        case 8: genotype[0]=2; genotype[1]=3; break;
//        case 9: genotype[0]=3; genotype[1]=3; break;
//        default: {printf("error 2 in genotype index conversion table (case=%i)",inde); exit(-1);}
//    }
//}

//populates a vector with the names of which we have data
//Should only open bcf file?
//Return data names
vector<char *> hasdata(char *fname){
    htsFile * inf = NULL;
    inf=hts_open(fname, "r");
    assert(inf);
    bcf_hdr_t *hdr = NULL;
    hdr=bcf_hdr_read(inf);
    assert(hdr);
    
    bcf1_t *rec = NULL;
    rec=bcf_init();
    assert(rec);
    hts_idx_t *idx=NULL;
    //idx=bcf_index_load(fname);
    idx=bcf_index_load(fname);
    assert(idx);
    hts_itr_t *iter=NULL;
    int nseq = 0;    // number of sequences
    const char **seqnames = NULL;
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    assert(seqnames);// bcf_hdr_id2name(hdr,i)
    
    vector<char *> ret;
    for(int i=0;i<nseq;i++){
        char buf[strlen(seqnames[i])+100];
        snprintf(buf,strlen(seqnames[i])+100,"%s:1-1000000000",seqnames[i]);
        iter=bcf_itr_querys(idx,hdr,buf);
        if(bcf_itr_next(inf, iter, rec)==0)
            ret.push_back(strdup(seqnames[i]));
        hts_itr_destroy(iter);
    }
    free(seqnames);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(inf);
    fprintf(stderr,"\t-> Done with preliminary parsing of file: we have data for %lu out of %d reference sequences\n",ret.size(),nseq);
    return ret;
}

const int PHREDMAX=256;
float pl2ln[PHREDMAX];

float pl2ln_f(int32_t& val){
    if(val>=PHREDMAX){
        return log(pow(10.0,-0.1*val));
    }else{
        return pl2ln[val];
    }
}

template <class T>
bool same_val_vcf(T a, T b) {
    return fabs(a - b) < numeric_limits<T>::epsilon();
}

bool is_nan_vcf(double x){
    return x != x;
}

//Modified from angsd
double emFrequency(double *loglike, int numInds, int iter, double start, char *keep, int keepInd){
    if(keepInd == 0){
        return 0.0;
    }
    float W[10];
    float p[4] = {(float)start,(float)start,(float)start,(float)start};
    float temp_p[4] = {(float)start,(float)start,(float)start,(float)start};
    double accu = 0.00001;
    double accu2 = 0;
    
    int it = 0;
    for(it=0;it<iter;it++){
        float sum[4]={0,0,0,0};
        for(int i=0;i<numInds;i++){
            if(keep!=NULL && keep[i]==0){
                continue;
            }
            int genotype[2];
            float sum2W = 0;
            for (int j=0;j<10;j++){
                findgenotypes_from_index(j, genotype);
                if(genotype[0] == genotype[1]){
                    W[j]=exp(loglike[i*10+j])*pow(p[genotype[0]],2);
                }else{
                    W[j]=exp(loglike[i*10+j])*2*p[genotype[0]]*p[genotype[1]];
                }
                sum2W += W[j];
            }
            //            cout<<sum2W<<"\n";
            if (sum2W == 0){
                fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) W:(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",i,loglike[i*10],loglike[i*10+1],loglike[i*10+2],loglike[i*10+3],loglike[i*10+4],loglike[i*10+5],loglike[i*10+6],loglike[i*10+7],loglike[i*10+8],loglike[i*10+9],W[0],W[1],W[2],W[3],W[4],W[5],W[6],W[7],W[8],W[9]);
            }
            assert(sum2W!=0);
            sum2W = 2*sum2W;
            sum[0] += (2*W[0]+W[1]+W[2]+W[3])/sum2W;
            sum[1] += (2*W[4]+W[1]+W[5]+W[6])/sum2W;
            sum[2] += (2*W[7]+W[2]+W[5]+W[8])/sum2W;
            sum[3] += (2*W[9]+W[3]+W[6]+W[8])/sum2W;
            if(isnan(sum[0])||isnan(sum[1])||isnan(sum[2])||isnan(sum[3])){
                fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) W:(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) sum=(%f,%f,%f,%f)\n",i,loglike[i*10],loglike[i*10+1],loglike[i*10+2],loglike[i*10+3],loglike[i*10+4],loglike[i*10+5],loglike[i*10+6],loglike[i*10+7],loglike[i*10+8],loglike[i*10+9],W[0],W[1],W[2],W[3],W[4],W[5],W[6],W[7],W[8],W[9],sum[0],sum[1],sum[2],sum[3]);
            }
        }
        float d = 0;
        for(int i=0;i<4;i++){
            p[i]=sum[i]/keepInd;
            d += pow(p[i]-temp_p[i],2);
        }
        if(sqrt(d)<accu){
            break;
        }
        for(int i=0;i<4;i++){
            temp_p[i]=p[i];
        }
    }
      if(isnan(p[0])||isnan(p[1])||isnan(p[2])||isnan(p[3])){
            return -999;
     }
    
    float k;
    for(int i=0;i<3;i++){
        for(int j=0;j<3-i;j++){
            if (temp_p[j]>temp_p[j+1]){
                k=temp_p[j];
                temp_p[j]=temp_p[j+1];
                temp_p[j+1]=k;
            }
        }
    }
    //Return the second latgest frequency
    return temp_p[2];
}

////from angsd
//double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){
//      if(keepInd == 0)
//            return 0.0;
//
//      float W0;
//      float W1;
//      float W2;
//      // fprintf(stderr,"start=%f\n",start);
//      float p=(float)start;
//      float temp_p=(float)start;
//      double accu=0.00001;
//      double accu2=0;
//      float sum;
//      int it=0;
//
//      for(it=0;it<iter;it++){
//            sum=0;
//            for(int i=0;i<numInds;i++){
//                  if(keep!=NULL && keep[i]==0)
//                        continue;
//                  W0=exp(loglike[i*3+0])*pow(1-p,2);
//                  W1=exp(loglike[i*3+1])*2*p*(1-p);
//                  W2=exp(loglike[i*3+2])*(pow(p,2));
//                  sum+=(W1+2*W2)/(2*(W0+W1+W2));
//                  if(isnan(sum))
//                fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f) W(%f,%f,%f) sum=%f\n",i,loglike[i*3],loglike[i*3+1],loglike[i*3+2],W0,W1,W2,sum);
//                }
//
//            p=sum/keepInd;
//            if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
//                  break;
//            temp_p=p;
//          }
//
//      if(isnan(p)){
//            fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
//            fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
//            //print_array(stderr,loglike,3*numInds);
//            fprintf(stderr,"keepList (nInd)\n");
//            //print_array(stderr,keep,numInds);
//            fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);
//
//            for(int ii=0;1&&ii<numInds;ii++){
//                  if(keep!=NULL && keep[ii]==1){
//                // fprintf(stderr,"1\t");
//                for(int gg=0;gg<3;gg++)
//                  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
//                fprintf(stderr,"\n");
//                      }
//                }
//            sum=0;
//            for(int i=0;i<numInds;i++){
//                  if(keep!=NULL && keep[i]==0)
//                        continue;
//                  W0=exp(loglike[i*3+0])*pow(1-p,2);
//                  W1=exp(loglike[i*3+1])*2*p*(1-p);
//                  W2=exp(loglike[i*3+2])*(pow(p,2));
//                  sum+=(W1+2*W2)/(2*(W0+W1+W2));
//                  fprintf(stderr,"[%s.%s():%d] p=%f W %f\t%f\t%f sum=%f loglike: %f\n",__FILE__,__FUNCTION__,__LINE__,p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
//                  break;
//                }
//            p=-999;
//            assert(p!=999);
//            return p;
//          }
//    
//      return(p);
//}

size_t getgls(char*fname,vector<double *> &mygl, vector<double> &freqs,int minind,double minfreq, string &vcf_format_field, string &vcf_allele_field,char *seek){
    
      for(int i=0;i<PHREDMAX;i++){
            pl2ln[i] = log(pow(10.0,-0.1*i));
          }
      htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);
      bcf_hdr_t *hdr = bcf_hdr_read(inf);
      bcf1_t *rec = NULL;rec=bcf_init();assert(rec);
      hts_idx_t *idx=NULL;
      hts_itr_t *iter=NULL;
    
      if(seek){
            fprintf(stderr,"\t-> Setting iterator to: %s\n",seek);fflush(stderr);
            idx=bcf_index_load(fname);
            iter=bcf_itr_querys(idx,hdr,seek);
          }
      
      int n    = 0;  // total number of records in file
      int nsnp = 0;  // number of SNP records in file
      int nseq = 0;  // number of sequences
      int nsamples = 0;
     // pl data for each call
      int npl_arr = 0;
      int npl     = 0;
      int32_t *pl = NULL;
     // gt data for each call
      int32_t ngt_arr = 0;
      int32_t ngt     = 0;
      int32_t *gt     = NULL;
    
    // af1/af data for each call
      int naf_arr = 0;
      int naf     = 0;
      float *af     = NULL;
    
    // read header
      nsamples = bcf_hdr_nsamples(hdr);
    
      const char **seqnames = NULL;
      seqnames = bcf_hdr_seqnames(hdr, &nseq); assert(seqnames);//bcf_hdr_id2name(hdr,i)
    
      char *chr;
      while(1){
            if(seek==NULL){
                  if(bcf_read(inf,hdr,rec)!=0)
                break;
                }else{
                      if(bcf_itr_next(inf, iter, rec)!=0)
                    break;
                    }
            n++;
            if (!bcf_is_snp(rec))
                  continue;
        nsnp++;
        if(rec->n_allele==1){
            continue;
        }
        int k1 = rec->n_allele*(rec->n_allele+1)/2;
        float ln_gl1[k1*nsamples];
        float ln_gl[10*nsamples];
        if(vcf_format_field == "PL") {
            for (int i = 0; i<10*nsamples;i++){
                ln_gl[i] = -INFINITY;
            }
                  npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
            
                  if(npl<0){
                        fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching PL tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
                        continue;
                      }
                  for (int i=0; i<npl; i++){
                        if ( pl[i]==bcf_int32_missing ){
                              bcf_float_set_missing(ln_gl1[i]);
                            } else if ( pl[i]==bcf_int32_vector_end ){
                                  bcf_float_set_vector_end(ln_gl1[i]);
                                } else{
                                      ln_gl1[i] = pl2ln_f(pl[i]);
                                    }
                      }
            
            int flag = 0;
            int nuc_num = 0;
            vector<int > renuc;
            vector<int > nuc;
            for (int n = 0; n < rec->n_allele; n++){
                int i = findnuctypeindex(rec->d.allele[n]);
                if (i == 4){
                    flag = 1;
                }else{
                    nuc.push_back(i);
                    nuc_num = nuc_num + 1;
                }
            }
            int renuc_num = 4 - nuc_num;
            if (renuc_num > 0 && flag == 1){
                for (int j = 0; j < 4; j++){
                    renuc.push_back(j);
                }
                for (int i = 0; i<nuc_num-1;i++){
                    for (int j = 0; j<nuc_num-1-i; j++){
                        if (nuc[j] > nuc[j+1]){
                            int k = nuc[j];
                            nuc[j] = nuc[j+1];
                            nuc[j+1] = k;
                        }
                    }
                }
                for (int i = nuc_num - 1; i > -1; i--){
                    renuc.erase(renuc.begin()+nuc[i]);
                }
            }
            for (int n1 = 0; n1 < rec->n_allele; n1++){
                int k = n1*(n1+1)/2;
                int i = findnuctypeindex(rec->d.allele[n1]);
                for (int n2 = 0; n2 <= n1; n2++){
                    int j = findnuctypeindex(rec->d.allele[n2]);
                    int gind = findgenotypeindex(i,j);
                    if (gind < 10){
                        for (int ns = 0; ns < nsamples; ns ++){
                            ln_gl[gind+ns*10] = ln_gl1[k+n2+ns*k1];
                        }
                    }else if (gind < 14){
                        for (int i1 = 0; i1 < renuc.size(); i1++){
                            int gind1 = findgenotypeindex(gind-10,renuc[i1]);
                            for (int ns = 0; ns < nsamples; ns ++){
                                ln_gl[gind1+ns*10] = ln_gl1[k+n2+ns*k1];
                            }
                        }
                    }else if (gind == 14){
                        for (int i1 = 0; i1 < renuc.size(); i1++){
                            int gind1 = findgenotypeindex(renuc[i1],renuc[i1]);
                            for (int ns = 0; ns < nsamples; ns ++){
                                ln_gl[gind1+ns*10] = ln_gl1[k+n2+ns*k1];
                            }
                        }
                    }
                }
            }
            for(int ns = 0; ns < nsamples; ns ++){
                float d = 0;
                for (int i = 0;i<k1;i++){
                    d = d + ln_gl1[i+ns*k1];
                }
                if (d == 0) {
                    for (int i = 0;i<10;i++){
                        ln_gl[i+ns*10]=0;
                    }
                }
            }
            
                } else if(vcf_format_field == "GT"){
                       int ngts = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
                       if ( ngts<0 ){
                             fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching GT tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
                             continue;
                           }
                       for(int ns=0; ns<nsamples;ns++){
                             int32_t *curr_ptr = gt + ns*2;
                    float *ln_gl_ptr = ln_gl + ns*10;
                             if ( bcf_gt_is_missing(curr_ptr[0]) || bcf_gt_is_missing(curr_ptr[1]) ){// obs genotype is missing
                        for (int i = 0;i<10;i++){
                            ln_gl_ptr[i] = 0;
                        }
                                 }else{
                            
                            for (int i = 0;i<10;i++){
                                ln_gl_ptr[i] = -INFINITY;
                            }
                            int al1 = bcf_gt_allele(curr_ptr[0]);
                            int al2 = bcf_gt_allele(curr_ptr[1]);
                            ln_gl_ptr[findgenotypeindex(al1,al2)]=0;
                        }
                         }
                
                   } else {
                          fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, GT or PL\n");
                       }
        int keepInd=0;
            char keep[nsamples];
            double *tmp = new double[10*nsamples];
            for(int ns=0;ns<nsamples;ns++){
                  float *ary= ln_gl+ns*10;
                  if ((is_nan_vcf(ary[0]) || is_nan_vcf(ary[1]) || is_nan_vcf(ary[2]) || is_nan_vcf(ary[3]) || is_nan_vcf(ary[4]) || is_nan_vcf(ary[5]) || is_nan_vcf(ary[6]) || is_nan_vcf(ary[7]) || is_nan_vcf(ary[8]) || is_nan_vcf(ary[9]) || ((ary[0] == 0) && (ary[1] == 0) && (ary[2] == 0) && (ary[3] == 0) && (ary[4] == 0) && (ary[5] == 0) && (ary[6] == 0) && (ary[7] == 0) && (ary[8] == 0) && (ary[9] == 0)))){
                        keep[ns]=0;
                      }else{
                    keep[ns]=1;
                    keepInd++;
                          }
            for(int i = 0;i<10;i++){
                      tmp[ns*10+i] = ary[i];
            }
                }
        
            
            naf = bcf_get_info_float(hdr, rec, vcf_allele_field.c_str(), &af, &naf_arr);
        
            double freq;
            if(naf==1){
                  freq = af[0];
                }else{
                      freq = emFrequency(tmp,nsamples,50,0.25,keep,keepInd);
                    }
            if(freq>0.999){
            freq=1;
        }else if(freq<0.001){
            freq=0;
        }
        if(keepInd>minind && freq>=minfreq && freq<= (1-minfreq)){
                  for(int ns=0;ns<10*nsamples;ns++){
                  tmp[ns]=exp(tmp[ns]);}
                  mygl.push_back(tmp);
                  freqs.push_back(freq);
                  //populate debug names
            
                } else {
                      delete [] tmp;
                    }
          }
      fprintf(stderr, "\t-> [file=\'%s\'][chr=\'%s\'] Read %i records %i of which were SNPs. Number of sites used for downstream analysis (MAF >= %f):%lu\n",fname,seek, n, nsnp, minfreq,mygl.size());
    
      free(pl);
      free(gt);
      bcf_hdr_destroy(hdr);
      bcf_close(inf);
      bcf_destroy(rec);
    
      if(iter)
            hts_itr_destroy(iter);
      if(idx)
            hts_idx_destroy(idx);
    
      //for(int i=0;i<nseq;i++)
      free(seqnames);
      return nsamples;
}

void *wrap(void *ptr){
      satan *god = (satan*) ptr;
      god->nind=getgls(god->fname, god->mygl,god->freqs, god->minind, god->minfreq,god->vcf_format_field,god->vcf_allele_field,god->seek);
      pthread_exit(NULL);//this is sometimes called without thread
}

int wrap_nothreading(void *ptr){
      satan *god = (satan*) ptr;
      god->nind=getgls(god->fname, god->mygl,god->freqs, god->minind, god->minfreq,god->vcf_format_field,god->vcf_allele_field,god->seek);
      return 1;
}

void *wrap2(void *){
      while(1){
            pthread_mutex_lock(&mymut);
            int myvar = mycounter;
            mycounter++;
            pthread_mutex_unlock(&mymut);
            if(myvar>=jobs.size())
                  pthread_exit(NULL);
            satan *god = (satan*) &jobs[myvar];
            god->nind=getgls(god->fname, god->mygl,god->freqs, god->minind, god->minfreq,god->vcf_format_field,god->vcf_allele_field,god->seek);
          }
}

double ** readbcfvcf(char*fname,int &nind, int &nsites, vector<double> &freqs,int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek){
      fprintf(stderr,"\t-> readbcfvcf seek:%s nind:%d\n",seek,nind);
      htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);
      bcf_hdr_t *hdr = NULL;hdr=bcf_hdr_read(inf);assert(hdr);
      int isbcf=0;
      vector<char *> hd;
    
      satan god;
      god.fname=fname;
      god.minind=minind;
      god.minfreq=minfreq;
      god.vcf_format_field=vcf_format_field;
      god.vcf_allele_field=vcf_allele_field;
      god.seek=seek;
    
      if(inf->format.format==bcf){
            isbcf=1;
            hd = hasdata(fname);
          }
    //hd basically tells you how many chromosomes are in the data.
      if(seek&&isbcf==0){
            fprintf(stderr,"\t-> if choosing region then input file has to be bcf\n");
            exit(0);
          }
    
      int nseq = 0;  // number of sequences
      const char **seqnames = NULL;
      seqnames = bcf_hdr_seqnames(hdr, &nseq);
      assert(seqnames);
      
      double **gls=NULL;
      if(seek!=NULL||isbcf==0){//single run
        
        wrap_nothreading(&god);
            nind=god.nind;
        
            gls=new double *[god.mygl.size()];
            for(int i=0;i<god.mygl.size();i++){
                  gls[i] = god.mygl[i];
                }
            freqs=god.freqs;
            nsites = god.mygl.size();
        fprintf(stderr,"Done reading everything in [%s] we have nsites:%d for samples:%d\n",seek,nsites,nind);
    }else{
            for(int i=0;i<hd.size();i++){
                  jobs.push_back(god);
                  jobs[i].seek=hd[i];
        }
            
            if(diskio_threads==1||isbcf==0){
                  for(int i=0;i<jobs.size();i++){
                wrap_nothreading(&jobs[i]);
                      }
                }else{
                      int at =0;
                      while(at<hd.size()){
                    int howmany=min(diskio_threads,(int)hd.size()-at);
                    pthread_t mythd[howmany];
                    for(int i=0;i<howmany;i++){
                                      if(std_queue)
                                              assert (pthread_create(&mythd[i],NULL,wrap,&jobs[i+at])==0);
                                     else
                                              assert (pthread_create(&mythd[i],NULL,wrap2,NULL)==0);
                    }
                    for(int i=0;i<howmany;i++)
                                assert(pthread_join(mythd[i],NULL)==0);
                    at+=howmany;
                          }
                    }
            nsites =0;
            for(int i=0;i<jobs.size();i++)
              nsites += jobs[i].mygl.size();
            nind=jobs[0].nind;
             fprintf(stderr,"Done reading everything we have nsites:%d for samples:%d\n",nsites,nind);
            //merge results
            gls=new double *[nsites];
            freqs.reserve(nsites);
            int at =0;
            for(int i=0;i<jobs.size();i++){
                  for(int j=0;j<jobs[i].mygl.size();j++)
            gls[at++] = jobs[i].mygl[j];
                  freqs.insert(freqs.end(),jobs[i].freqs.begin(),jobs[i].freqs.end());
                }
            for(int i=0;i<hd.size();i++)
              free(hd[i]);
          }
    free(seqnames);
      if(hdr) bcf_hdr_destroy(hdr);
      hts_close(inf);
      return gls;
}

double vcftwoDSFS(char* fname, int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek){
    int nsites = 0;
    double **GLDATA=NULL;
    int nind;
    vector<double> freqs;
    double t, parameters[8], twoDSFS[10][10];
    GLDATA = readbcfvcf(fname,nind,nsites,freqs,minind,minfreq,vcf_format_field,vcf_allele_field,seek);
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    cout<<nsites<<"\n";
    
    int *SDATA = (int *) malloc(nsites * sizeof(int));
    
    int eff_nsites = FilterSites(GLDATA, SDATA, nsites);
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_nsites, nsites);
    
    //Estimate T
    estimateT(twoDSFS, &t, parameters);
    cout<<"Estimated t = "<<t<<"\n";
    
    for (int i = 0; i < nsites; i++){
        free(GLDATA[i]);
    }
    free(GLDATA);
    
    return t;
}

double vcfjointEM(char* fname, int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek, double* tt){
    int nsites = 0;
    double **GLDATA=NULL;
    int nind;
    vector<double> freqs;
    double t, parameters[8], twoDSFS[10][10];
    
    GLOBPsum = (double**) malloc(40 * sizeof(double ***));
    for (int i = 0; i < 40; i++){
        GLOBPsum[i] = (double *) malloc(40 * sizeof(double));
    }
    
    GLDATA = readbcfvcf(fname,nind,nsites,freqs,minind,minfreq,vcf_format_field,vcf_allele_field,seek);
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    int *SDATA = (int *) malloc(nsites * sizeof(int));
    
    int eff_numsites = FilterSites(GLDATA, SDATA, nsites);
    
    EMlikelihoodforT(GLDATA, SDATA, nsites, 0.001, 0.00004, 0.0002,&tt[0]);
    
    t = tt[0]+tt[1]+tt[2];
    cout<<"Estimated t = "<<t<<"\n";
    
//    for (int i = 0; i < nsites; i++){
//        free(GLDATA[i]);
//    }
//    free(GLDATA);
    
//    for (int i = 0; i < 40; i++){
//        free(GLOBPsum[i]);
//    }
//    free(GLOBPsum);
    
    return t;
}

//#ifdef __WITH_MAIN__

//int main(){
//    char* name = "lei.bcf";
//    int nsites = 0;
//    double **gls=NULL;
//    int nind;
//    vector<double> freqs;
//    string pl=string("PL");
//    string fr=string("AFngsrelate");
//    char *reg = NULL;
//    gls = readbcfvcf(name,nind,nsites,freqs,1,0.45,pl,fr,reg);
//    //The first index is how many site (seems to filter out multialleltic sites) across different chromosomes.
//    //The second index is 3*number of individuals.
//    //Both of these should be modified.
//    for (int i = 0;i<nsites; i++){
//        for (int j = 0;j<20;j++){
//            cout<<gls[i][j]<<" ";
//        }
//        cout<<"\n";
//    }
//    for (int i = 0;i<nsites; i++){
//        cout<<freqs[i]<<"\n";
//    }
//    //    cout<<freqs[0]<<" "<<freqs[1]<<"\n";
//    //    cout<<gls[0][0]<<" "<<gls[0][1]<<" "<<gls[0][2]<<" "<<gls[0][3]<<" "<<gls[0][4]<<" "<<gls[0][5]<<" "<<gls[0][6]<<" "<<gls[0][7]<<" "<<gls[0][8]<<"\n";
//    //    cout<<gls[1][0]<<" "<<gls[1][1]<<" "<<gls[1][2]<<" "<<gls[1][3]<<" "<<gls[1][4]<<" "<<gls[1][5]<<" "<<gls[1][6]<<" "<<gls[1][7]<<" "<<gls[1][8]<<"\n";
//    return 0;
//}
// #endif

