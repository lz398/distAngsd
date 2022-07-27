//
//  mpileuptest.hpp
//  
//
//  Created by Lei Zhao on 24/11/2021.
//

#ifndef mpileuptest_h
#define mpileuptest_h
#include <vector>
#include <string>
#include <htslib/kstring.h>
#include <zlib.h>

std::vector<std::string> split(char* c, const char* src);
void op(std::vector<std::string> &src, const char* opc);
char* join(std::vector<std::string> src);
char* readproc(char* c, int &idep, char*& qscores);
int mpileupreader(gzFile gz, std::vector<kstring_t> &indiv0, std::vector<kstring_t> &indiv1, std::vector<kstring_t> &bq0, std::vector<kstring_t> &bq1);
int mpileuptwoDSFS(gzFile gz, double par[9], int isthreading, int isuchar, int dobinary, int is2Dinfer, int isex, double &t, double &p);
void EMStepforNuc2DSFS_mpileup(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], std::vector<kstring_t> *indiv0, std::vector<kstring_t> *indiv1, std::vector<kstring_t> *bq0, std::vector<kstring_t> *bq1, size_t start, size_t numsites, size_t &psum);
double estimateNuc2DSFS_EM_mpileup(double SEQ2DSFS[4][4], std::vector<kstring_t> *indiv0, std::vector<kstring_t> *indiv1, std::vector<kstring_t> *bq0, std::vector<kstring_t> *bq1, size_t numsites);
#endif /* mpileuptest_hpp */
