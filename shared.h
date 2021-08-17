//
//  special.h
//  
//
//  Created by Lei Zhao on 14/01/2021.
//

#ifndef special_h
#define special_h
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
#include <htslib/vcf.h>
#include <limits>
#include <pthread.h>
#include <cassert>
using namespace Eigen;
using namespace std;



extern double tole;

extern MatrixXd RRVEC, LRVEC, PMAT, PMATdt, PMATdt2;
extern VectorXd RRVAL, pi; /*RRVEC and LRVEC are globals used by
                           diagonalizeGTR and gettransitionprobmatGTR.
                           PMAT and pi is the transition probability matrix and stationary distcb., respectively, accessed by diagonalizeGTR and gettransitionprobmatGTR and pijtGTR.
                           */
extern double GLOBtwoDSFS[10][10], GLOBSEQ2DSFS[4][4], GLOBpar[8], globerror, globt1, globt2;
extern size_t globnumsites;
extern double **SEQDATA, **GLOBPsum;
extern int *SEQ_SDATA;
extern int **DATA; /*the data is stored here*/
extern int newt, newmat; /*these are control variables used to determine if the transition probability matrix needs to be recalculated*/
const int PHREDMAX=256;
extern float pl2ln[PHREDMAX];
extern double plmatrix[PHREDMAX][PHREDMAX];
typedef unsigned char uchar;
struct EMjob{
    size_t index;
    double **gls;
    int *ss;
    size_t start;
    size_t len;
    int rowlen;
    double segESFS2[10][10];
    double two2DSFS[10][10];
};
struct EMjob_uchar{
    size_t index;
    uchar **gls;
    size_t start;
    size_t len;
    int rowlen;
    double segESFS2[10][10];
    double two2DSFS[10][10];
};
struct double4{
    double vec[4];
};
struct EMjobforSEQ2DSFS{
    size_t index;
    size_t psum;
    vector<vector<double4> > p0;
    vector<vector<double4> > p1;
    size_t start;
    size_t len;
    double segESEQSFS2[4][4];
    double SEQ2DSFS[4][4];
};
#endif /* special_h */
