#ifndef GLtest_h
#define GLtest_h

#include<htslib/bgzf.h>
#include<htslib/kstring.h>
#include "shared.h"

using namespace std;
using namespace Eigen;

/*Randomness Generators*/
int static z_rndu=137;
void SeedSetup();
void SetSeed(int seed);
// gamma distribution
double gammln(double xx);
// Poisson distribution
double Poisson(double xm);
double uniform();
double uniform_yang();

/*Different substitution models*/
//hard coded transition probabilities under JC69 model.  '*par' is not used as JC models uses no parameters beyond t.
double pijtJC(double t, double *par, int from, int to);
//hard coded transition probabilities under HKY85 model
//This code has not been tested
double pijtHKY(double t, double *par, int from, int to);
//hard coded transition probabilities under F84 model
//this code has not been tested
double pijtF84(double t, double *par, int from, int to);
void diagonalizeGTR( double *par);
void gettransitionprobmatGTR(double t);
double pijtGTR(double t, double *par,  int from, int to);
//double like(int numsites, int **data, double (*pijt)(double t, double *par, int from, int to), double *locpi, double t, double *par);


//Might be put into io.h
int readseqnuc();

/* General simulation*/
void simnucleotides(int Genotype[2], double simmat[4][4]);
void makeSIMMAT(double simmat[4][4]);
void makepolyMAT(double simmat[4][4]);
void simpoly(int nuc, double SIMMAT[4][4], int genotype[2]);
int findgenotypeindex(int i, int j);
void findgenotypes_from_index(int inde, int genotype[2]);

//Not used
void makeGLswitchmatrix(double e, double switchmatrix[10][10], double simswitchmatrix[10][10]);
int simSEQs_per_read(int genotypes[2], double e);
void simSEQs_reads(int **SEQDATA, int genotypes[2], size_t site, int species, double e, double RD);


/* 1.2. Joint Genotype distribution method: standard and tree structure*/
//Simulation: Generate a read of a true genotype with calling error rate e, and update the genotype likehood.
int simGLs_per_read(double **GLDATA, int genotypes[2], size_t site, int species, double e);
/* Simulation: Generate Read depth for each site*/
int ReaddepthGenerator(double RD, int par);
// Simulation: Gather different reads and finalize the genotype likelihood for one site.
void simGLs_reads(double **GLDATA, int genotypes[2], size_t site, int species, double e, double RD);
//Simulation: Simulate two individuals' genome and reads in two species.
/*tdiv is the divergence time. t1 and t2 are the average coalescence times within species*/
void simulateGLsTwoSpecies(double RD, size_t numsites, double errorate, double tdiv, double t1, double t2, double **GLDATA,  double (*pijt)(double t, double *par, int from, int to), double *par);
//Filter: All sites considered should have data for both individual.
size_t FilterSites(double **GLDATA, int *SDATA, size_t numsites);
//Inference: The following two functions are used to measure the convergence of joint-distribution of genotypes.
void differr2D(double* m1, double* m2, double* mdiff, int n);
double err2D(double* m, int n);
//Inference: One EM step for joint genotype distribution
//void EMStepfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int numsites);
void EMStepfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int *SDATA, size_t numsites);
// Inference: Accelerate EM step for joint genotype distribution
//int EMAccelfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int numsites);
int EMAccelfor2DSFS(double twoDSFS[10][10], double ESFS2[10][10], double **GLDATA, int *SDATA, size_t numsites);
//Inference: Main EM algorithm for joint genotype distribution
//double estimate2DSFS_EM(double twoDSFS[10][10], double **GLDATA, int numsites);
double estimate2DSFS_EM(double twoDSFS[10][10], double **GLDATA, int *SDATA, size_t eff_numsites, size_t numsites);
double seqlike(int observedbase, int truebase, double errorrate);
//Inference: total likelihood based on the inferred joint genotype distribution.
double likeGLwithtwoDSFS(double twoDSFS[10][10], double t, double par[8]);
//Inference: total likelihood based on the inferred joint genotype distribution.
//Add in the tree structure
double likeGLwithtwoDSFS_m(double twoDSFS[10][10], double t, double par[8], double t1, double t2);
//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT(double t);
//Inference: Calculate the likelihood for divergence t, joint genotype distribution + tree structure
double likelihoodforT_m(double t);
//Inference: Estimation of divergence time t based on joint genotype distribution
double estimateT(double twoDSFS[10][10], double *t, double parameters[],int isex);
/*Inference: Estimation of divergence time t based on joint genotype distribution + tree structure.*/
double estimateT_m(double twoDSFS[10][10], double *t, double parameters[], double t1, double t2);
/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution*/
double testtwoDSFS(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double par[9], const char* glfname, int isthreading, int dobinary, int isex, int r);
/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + tree structure*/
double testtwoDSFS_m(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate);


/* 3. A Random Read method */
/*Simulation: Generate a read of a true genotype with calling error rate e, and based on the chosen read update likelihood.*/
int simSEQs_per_read_v1(double **SEQDATA, int genotypes[2], size_t site, int species, double e, int flag);
/*Simulation: Gather different reads and finalize the chosen read likelihood for one site. In the future should be able to switch from Poisson distributed read depth to constant coverage, and deal with 0 read depth.*/
void simSEQs_reads_v1(double **SEQDATA, int genotypes[2], size_t site, int species, double e, double RD);
/*Simulation: Simulate two individuals' genome and reads in two species.
 tdiv is the divergence time. t1 and t2 are the average coalescence times within species, based on the chosen read*/
void simulateGLsTwoSpeciesSEQ_v1(double RD, size_t numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par);
size_t FilterSitesSEQ(double **SEQDATA, int *SEQ_SDATA, size_t numsites);
//Inference: total likelihood based on a random chosen/consensus(future) read per site
//double SEQlikelihood_v1(double **SEQDATA, double t, double *par, int numsites);
double SEQlikelihood_v1(double **SEQDATA, int *SEQ_SDATA, double t, double *par, size_t numsites);
/*Inference: Calculate the likelihood based on a random chosen/consensus(future) read per site*/
double likelihoodforTSEQ_v1(double t);
/*Simulation+Inference: Simulation and estimation of divergence time t based on a chosen read per site*/
double testsimSEQDATA_v1(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate);


/* 4. Joint nucleotide distribution method */
/*A structure only for Joint nucleotide distribution method*/

/*Simulation: Generate a read of a true genotype with calling error rate e*/
int simSEQs_per_read_v2(int genotypes[2], size_t site, int species, double e);
/*Simulation: Gather different reads and calculate nucleotide likelihood per read for one site. In the future should be able to switch from Poisson distributed read depth to constant coverage, and deal with 0 read depth.*/
void simSEQs_reads_v2(vector<vector<double4> > &P, int genotypes[2], size_t site, int species, double e, double RD);
/*Simulation: Simulate two individuals' genome and reads in two species.
 tdiv is the divergence time. t1 and t2 are the average coalescence times within species, based on the joint selected nucleotide distribution*/
void simulateGLsTwoSpeciesSEQ_v2(double RD, size_t numsites, double errorrate, double tdiv, double t1, double t2, vector<vector<double4> >&P0, vector<vector<double4> >&P1, double (*pijt)(double t, double *par, int from, int to), double *par);
/*Check the effective number of sites*/
size_t CheckSites(vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites);
/*Inference: total likelihood based on the inferred joint selected nucleotide distribution.*/
double likeSEQwithtwoDSFS(double SEQ2DSFS[4][4], double t, double par[8]);
/*Inference: Calculate the likelihood for divergence t, joint selected nucleotide distribution*/
double likelihoodforTSEQ2DSFS(double t);
/*Inference: Estimation of divergence time t based on joint selected nucleotide distribution*/
double estimateTSEQ2DSFS(double SEQ2DSFS[4][4], double *t, double parameters[],int isex);
//Inference: One EM step for joint selected nucleotide distribution
void EMStepforNuc2DSFS(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites, size_t &psum);
// Inference: Accelerate EM step for joint selected nucleotide distribution
int EMAccelforNuc2DSFS(double SEQ2DSFS[4][4], double ESEQSFS2[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites);
//Inference: Main EM algorithm for joint genotype distribution
double estimateNuc2DSFS_EM(double SEQ2DSFS[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites);
/* Simulation + Inference: Simulation and estimation of divergence time t based on joint selected nucleotide distribution*/
double testsimSEQ2DSFS(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double par[9], int isthreading, int isex);


/* 5. Whole EM method for Inference*/
//void GetEMMatrix(double tdiv0, double t10, double t20, int numsites, double** Psum, double** GLDATA);
void GetEMMatrix(double tdiv0, double t10, double t20, size_t numsites, double** Psum, double** GLDATA, int *SDATA);
double Q1(double t1, int l);
static int globl1=1;
static int globl2=2;
double Q11(double t1);
double Q12(double t1);
double Q2(double tdiv);
//void EMStepforT(double* t0, double *t, int numsites,double **GLDATA);
void EMStepforT(double* t0, double *t, size_t numsites,double **GLDATA,int *SDATA);
void differr(double *v1, double *v2, double *v3, int dim);
double err(double *v, int dim);
//int EMAccelforT(double* t0, double *t, int numsites,double **GLDATA);
int EMAccelforT(double* t0, double *t, size_t numsites,double **GLDATA, int *SDATA);
//void EMlikelihoodforT(double **GLDATA, int numsites, double tdiv0, double t10, double t20, double* t);
void EMlikelihoodforT(double **GLDATA, int* SDATA, size_t numsites, double tdiv0, double t10, double t20, double* t);
double testjointEM(double RD, size_t numsites, double tdiv, double t1, double t2, double errorrate, double* tt);

/*6. threading EM for 2DSFS algorithm*/
void *athread(void *aptr);
void EMStepfor2DSFS_threading_initial(size_t numsites,int rowL,int nthreads,vector<EMjob> &jobvec,double **GLDATA, int *SDATA);
void EMStepfor2DSFS_threading(pthread_t* mythd, int nthreads,vector<EMjob> &jobvec,double **GLDATA,double two2DSFS[10][10], size_t numsites, size_t eff_numsites);
double estimate2DSFS_EM_threading(double twoDSFS[10][10], double **GLDATA, int* SDATA, size_t numsites, size_t eff_numsites, int nthreads, int rowL);

void *bthread(void *aptr);
void EMStepforNuc2DSFS_threading_initial(size_t numsites, int nthreads,vector<EMjobforSEQ2DSFS> &jobvec, vector<vector<double4> > &P0, vector<vector<double4> > &P1);
void EMStepforNuc2DSFS_threading(pthread_t* mythd, int nthreads,vector<EMjobforSEQ2DSFS> &jobvec,double SEQ2DSFS[4][4], size_t numsites);
double estimateNuc2DSFS_EM_threading(double SEQ2DSFS[4][4], vector<vector<double4> >&P0, vector<vector<double4> >&P1, size_t numsites,int nthreads);

int gls_writer_double(const char* glfname, int dobinary, int nsites, double **gls);
int gls_writer_uchar(const char* glfname, int dobinary, int nsites, uchar **gls);
#endif
