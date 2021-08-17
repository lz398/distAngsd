//
//  io.h
//  
//
//  Created by Lei Zhao on 12/08/2021.
//

#ifndef io_h
#define io_h

#include <stdio.h>
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
using namespace std;

typedef struct{
  //filenames
    const char* outname;
    const char* method;
    const char* model;
    const char* glfname;
    // Specify the GTR model
    double par[9];
    int dobinary;
    int is2Dinfer;
    int isthreading;
    
    // For simulation
    int simrep;
    double errorrate;
    int numsites;
    double RD;
    double tdiv;
    double t1;
    double t2;
    double p_inv;
    
    // For real data
    char* vcfname;
    int isuchar;
}pars;

pars *pars_init();

pars *get_pars(int argc,char **argv);

#endif /* io_h */
