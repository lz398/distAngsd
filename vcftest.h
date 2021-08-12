//
//  vcftest.h
//  
//
//  Created by Lei Zhao on 06/01/2021.
//

#ifndef vcftest_h
#define vcftest_h
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <cmath>
#include <limits>
#include <string>
#include <pthread.h>
#include <cassert>

using namespace std;

typedef struct satan_t{
    char*fname;
    int minind;
    double minfreq;
    string vcf_format_field;
    string vcf_allele_field;
    char *seek;
    vector<double *> mygl;
    vector<double> freqs;
    int nind;
}satan;
//vector<satan> jobs;

int findnuctypeindex(char* c1);
//int findgenotypeindex(int i, int j);
//void findgenotypes_from_index(int inde, int genotype[2]);
vector<char *> hasdata(char *fname);
float pl2ln_f(int32_t& val);
template <class T>
bool same_val_vcf(T a, T b);
bool is_nan_vcf(double x);
double emFrequency(double *loglike, int numInds, int iter, double start, char *keep, int keepInd);
size_t getgls(char*fname,vector<double *> &mygl, vector<double> &freqs,int minind,double minfreq, string &vcf_format_field, string &vcf_allele_field,char *seek);
void *wrap(void *ptr);
int wrap_nothreading(void *ptr);
void *wrap2(void *);
double ** readbcfvcf(char*fname,int &nind, int &nsites, vector<double> &freqs,int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek);

//vcf inference
double vcftwoDSFS(char* fname, int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek);
double vcfjointEM(char* fname, int minind,double minfreq, string vcf_format_field, string vcf_allele_field,char *seek, double* tt);
#endif /* vcftest_h */
