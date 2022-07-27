//
//  vcftest.h
//  
//
//  Created by Lei Zhao on 06/01/2021.
//

#ifndef vcftest_h
#define vcftest_h
#include <htslib/vcf.h>
#include "shared.h"

typedef struct satan_double_t{
    char*fname;
    int minind;
    double minfreq;
  std::string vcf_format_field;
  std::string vcf_allele_field;
    char *seek;
  std::vector<double *> mygl;
  std::vector<double> freqs;
    int nind;
}satan_double;

typedef struct satan_uchar_t{
    char*fname;
    int minind;
    double minfreq;
  std::string vcf_format_field;
  std::string vcf_allele_field;
    char *seek;
    //vector<double *> mygl;
    /*This is to save some precious space! The reason that this can be done since the genotype
     likelihoods in vcf files are of low resolution anyway.*/
  std::vector<uchar *> mygl;
  std::vector<double> freqs;
    int nind;
}satan_uchar;

//vector<satan> jobs;

int findnuctypeindex(char* c1);
std::vector<char *> hasdata(char *fname);
float pl2ln_f(int32_t& val);
template <class T>
bool same_val_vcf(T a, T b);
bool is_nan_vcf(double x);

double emFrequency_double(double *loglike, int numInds, int iter, double start, char *keep, int keepInd);
double emFrequency_uchar(uchar *pl_gl, int numInds, int iter, double start, char *keep, int keepInd);
size_t getgls_double(char*fname,std::vector<double *> &mygl, std::vector<double> &freqs,int minind,double minfreq, std::string &vcf_format_field, std::string &vcf_allele_field,char *seek);
size_t getgls_uchar(char*fname,std::vector<uchar *> &mygl, std::vector<double> &freqs,int minind,double minfreq, std::string &vcf_format_field, std::string &vcf_allele_field,char *seek);
void *wrap_double(void *ptr);
void *wrap_uchar(void *ptr);
int wrap_nothreading_double(void *ptr);
int wrap_nothreading_uchar(void *ptr);
void *wrap2_double(void *);
void *wrap2_uchar(void *);
double ** readbcfvcf_double(char*fname,int &nind, size_t &nsites, std::vector<double> &freqs,int minind,double minfreq, std::string vcf_format_field, std::string vcf_allele_field,char *seek);
uchar ** readbcfvcf_uchar(char*fname,int &nind, size_t &nsites, std::vector<double> &freqs,int minind,double minfreq, std::string vcf_format_field, std::string vcf_allele_field,char *seek);

/*7. threading EM for 2DSFS algorithm **double and **uchar genotype likelood*/
void EMStepfor2DSFS_uchar(double twoDSFS[10][10], double ESFS2[10][10], uchar **GLDATA, size_t numsites);
double estimate2DSFS_EM_uchar(double twoDSFS[10][10], uchar **GLDATA, size_t numsites);
void *athread_uchar(void *aptr);
void EMStepfor2DSFS_threading_initial_uchar(size_t numsites,int rowL,int nthreads,std::vector<EMjob_uchar> &jobvec,uchar **GLDATA);
void EMStepfor2DSFS_threading_uchar(pthread_t* mythd, int nthreads,std::vector<EMjob_uchar> &jobvec,uchar **GLDATA,double two2DSFS[10][10], size_t numsites);
double estimate2DSFS_EM_threading_uchar(double twoDSFS[10][10], uchar **GLDATA, size_t numsites, int nthreads, int rowL);

void EMStepfor2DSFS_double(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, size_t numsites);
double estimate2DSFS_EM_double(double twoDSFS[10][10], double **GLDATA, size_t numsites);
void *athread_double(void *aptr);
void EMStepfor2DSFS_threading_initial_double(size_t numsites,int rowL,int nthreads,std::vector<EMjob> &jobvec,double **GLDATA);
void EMStepfor2DSFS_threading_double(pthread_t* mythd, int nthreads,std::vector<EMjob> &jobvec,double **GLDATA,double two2DSFS[10][10], size_t numsites);
double estimate2DSFS_EM_threading_double(double twoDSFS[10][10], double **GLDATA, size_t numsites, int nthreads, int rowL);

//vcf inference
//double vcftwoDSFS(char* fname, const char* glfname, int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek, double par[9], int isthreading, int isuchar, int dobinary);
void vcftwoDSFS(char* fname, const char* glfname, int minind,double minfreq, std::string vcf_format_field, std::string vcf_allele_field, char *seek, double par[9], int isthreading, int isuchar,int dobinary, int is2Dinfer, int isex, double &t, double &p);

int CheckTable(const char* tabname, size_t &numsites, size_t &cols);
int CheckTablebin(const char* tabname, size_t &numsites, size_t &cols, size_t ind, int isuchar);
void gls_read(const char* tabname, double **GLDATA, int isuchar, int dobinary);
void tabletwoDSFS(const char* tabname, double par[9], int isthreading, int inuchar, int inbinary, int is2Dinfer, int isex, double &t, double &p);
#endif /* vcftest_h */
