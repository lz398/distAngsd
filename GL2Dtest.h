//
//  GL2Dtest.h
//  
//
//  Created by Lei Zhao on 01/06/2021.
//

#ifndef GL2Dtest_h
#define GL2Dtest_h
#include "shared.h"

void simnucleotidesinv(int Genotype[2],double simmat[4][4]);
void simulateGLsTwoSpeciesWithInvSite(double RD, size_t numsites, double p_inv, double errorate, double tdiv, double t1, double t2, double **GLDATA,  double (*pijt)(double t, double *par, int from, int to), double *par);
void testtwoDSFSWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9], const char* glfname, int isthreading, int dobinary, int r);
double likeGLwithtwoDSFSWithInvSite(double twoDSFS[10][10], const double* x, double par[8]);
double likelihoodforTWithInvSite(const double* x,const void *);
double estimateTWithInvSite(double twoDSFS[10][10], double x[2], double parameters[]);
void likeGLwithtwoDSFSWithInvSite_gradient(double twoDSFS[10][10], const double* x, double *y, double par[8]);
void likelihoodforTWithInvSite_gradient(const double* x, double* y);
void gettransitionprobmatGTR_gradient(double t);
void simulateGLsTwoSpeciesSEQ2DSFSWithInvSite(double RD, size_t numsites, double p_inv, double errorrate, double tdiv, double t1, double t2, vector<vector<double4> >&P0, vector<vector<double4> >&P1, double (*pijt)(double t, double *par, int from, int to), double *par);
double likeSEQwithtwoDSFSWithInvSite(double SEQ2DSFS[4][4], const double* x, double par[8]);
double likelihoodforTSEQ2DSFSWithInvSite(const double* x, const void *);
void likeSEQwithtwoDSFSWithInvSite_grad(double SEQ2DSFS[4][4], const double* x, double *y, double par[8]);
void likelihoodforTSEQ2DSFSWithInvSite_grad(const double* x, double* y);
double estimateTSEQ2DSFSWithInvSite(double SEQ2DSFS[4][4], double x[2], double parameters[]);
void testsimSEQ2DSFSWithInvSite(double RD, int numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9], int isthreading);
void testsimSEQDATA_randomWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9]);
void testsimSEQDATA_consensusWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p, double par[9]);
void testGL_consensusWithInvSite(double RD, size_t numsites, double p_inv, double tdiv, double t1, double t2, double errorrate, double &t, double &p);
void findmaxvec(double vec[], int n, vector<int> &indices);

double estimateTWithInvSitePinvKnown(double twoDSFS[10][10], double *t, double parameters[],double p_inv);
double likelihoodforTWithInvSitePinvKnown(double t);
double estimateTSEQ2DSFSWithInvSitePinvKnown(double SEQ2DSFS[4][4], double* t, double parameters[],double p_inv);
double likelihoodforTSEQ2DSFSWithInvSitePinvKnown(double t);
//
double parlikeGLwithtwoDSFSWithInvSite(double twoDSFS[10][10], double t, double par[8]);
double parlikelihoodforTWithInvSite(double t);
double estimateTWithInvSiteParlike(double twoDSFS[10][10], double *t, double parameters[]);
double parlikeSEQwithtwoDSFSWithInvSite(double SEQ2DSFS[4][4], const double t, double par[8]);
double parlikelihoodforTSEQ2DSFSWithInvSite(double t);
double estimateTSEQ2DSFSWithInvSiteParlike(double SEQ2DSFS[4][4], double* t, double parameters[]);
#endif /* GL2Dtest_h */
