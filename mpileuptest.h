//
//  mpileuptest.hpp
//  
//
//  Created by Lei Zhao on 24/11/2021.
//

#ifndef mpileuptest_h
#define mpileuptest_h

#include <stdio.h>
#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <htslib/kstring.h>
using namespace std;
#define LENS 20480
vector<string> split(char* c, const char* src);
void op(vector<string> &src, const char* opc);
char* join(vector<string> src);
char* readproc(char* c, int &idep, char*& qscores);
int mpileupreader(gzFile gz, vector<kstring_t> &indiv0, vector<kstring_t> &indiv1, vector<kstring_t> &bq0, vector<kstring_t> &bq1);
int mpileuptwoDSFS(gzFile gz, double par[9], int isthreading, int isuchar, int dobinary, int is2Dinfer, double &t, double &p);
void EMStepforNuc2DSFS_mpileup(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<kstring_t> *indiv0, vector<kstring_t> *indiv1, vector<kstring_t> *bq0, vector<kstring_t> *bq1, size_t start, size_t numsites, size_t &psum);
double estimateNuc2DSFS_EM_mpileup(double SEQ2DSFS[4][4], vector<kstring_t> *indiv0, vector<kstring_t> *indiv1, vector<kstring_t> *bq0, vector<kstring_t> *bq1, size_t numsites);
#endif /* mpileuptest_hpp */
