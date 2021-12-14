//
//  NewMethods.hpp
//  
//
//  Created by Lei Zhao on 15/01/2021.
//

#ifndef NewMethods_hpp
#define NewMethods_hpp
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
#include "GLtest.h"

//Inference: total likelihood based on the inferred joint genotype distribution + _Ambiguity
double likeGLwithtwoDSFS_Ambiguity(double twoDSFS[10][10], double t, double par[8]);
//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT_Ambiguity(double t);
//Inference: Estimation of divergence time t based on joint genotype distribution + _Ambiguity
double estimateT_Ambiguity(double twoDSFS[10][10], double *t, double parameters[]);
/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + _Ambiguity*/
double testtwoDSFS_Ambiguity(double RD, int numsites, double tdiv, double t1, double t2, double errorrate);

//Inference: total likelihood based on the inferred joint genotype distribution + _Ambiguity
double likeGLwithtwoDSFS_Ambiguity_v2(double twoDSFS[10][10], double t, double par[8]);
//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT_Ambiguity_v2(double t);
//Inference: Estimation of divergence time t based on joint genotype distribution + _Ambiguity
double estimateT_Ambiguity_v2(double twoDSFS[10][10], double *t, double parameters[]);
/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + _Ambiguity*/
double testtwoDSFS_Ambiguity_v2(double RD, int numsites, double tdiv, double t1, double t2, double errorrate);
//double gen_dist(double twoDSFS[10][10]);
//void gen_dist(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites);
void gen_dist(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites, double &dist, double &dist_JC, double &dist1, double &dist1_JC);
double testtwoDSFS_ngsDist(double RD, int numsites, double tdiv, double t1, double t2, double errorrate);
// Another Extension of ngsDist: Basically only account for "SNP"
void gen_dist_snp(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites, double &dist, double &dist_JC, double &dist1, double &dist1_JC);
double testtwoDSFS_ngsDist_snp(double RD, int numsites, double tdiv, double t1, double t2, double errorrate);

void simSEQs_reads_consensus(double **SEQDATA, int genotypes[2], int site, int species, double e, double RD);
void simSEQs_reads_random(double **SEQDATA, int genotypes[2], int site, int species, double e, double RD);
void estimate2DSFS_consensusGT(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites);
void simulateGLsTwoSpeciesSEQ_random(double RD, int numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par);
void simulateGLsTwoSpeciesSEQ_consensus(double RD, int numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par);
// RandomSEQ
double testsimSEQDATA_random(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9]);
// ConsensusSEQ
double testsimSEQDATA_consensus(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9]);
// ConsensusGT
double testtwoDSFS_consensusGT(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9]);
#endif /* NewMethods_hpp */
