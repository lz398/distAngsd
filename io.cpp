//
//  io.cpp
//
//
//  Created by Lei Zhao on 12/08/2021.
//

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
#include "io.h"

using namespace std;

pars *pars_init(){
    pars *p =(pars*) calloc(1,sizeof(pars));
    p->outname = strdup("distAngsdlog");
    // method can be geno, nuc, RandomSEQ, ConsensusSEQ, ConsensusGT for simulation
    // method can only be geno for vcf
    p->method = strdup("geno");
    // model can be either JC or GTR
    p->model = strdup("JC");
    p->glfname = NULL;
    p->vcfname = NULL;
    p->dobinary = 1;
    p->simrep = NULL;
    p->is2Dsim = 0;
    p->p_inv = 0;
    for (int i=0;i<5;i++){
        p->par[i] = 1.0;
    }
    for (int i=5;i<9;i++){
        p->par[i] = 0.25;
    }
    
    p->errorrate = 0.002;
    p->numsites = 1000000;
    p->RD = 1;
    p->tdiv = 1;
    p->t1 = 0.4;
    p->t2 = 0.25;
    p->isthreading = 0;
    p->isuchar = 0;
    return p;
}

pars *get_pars(int argc,char **argv){
    pars *p = pars_init();
    if(argc % 2){
        fprintf(stderr,"\t-> Must supply arguments in the form -pattern value\n");
        free(p);
        return NULL;
    }
    int k = 0;
    double sumpi = 0;
    while(*argv){
        char *key=*argv;
        char *val=*(++argv);
        if(!strcasecmp("-o",key)) p->outname=strdup(val);
        else if(!strcasecmp("-method",key)) p->method=strdup(val);
        else if(!strcasecmp("-model",key)) p->model=strdup(val);
        else if(!strcasecmp("-glf",key)) p->glfname=strdup(val);
        else if(!strcasecmp("-vcf",key)) p->vcfname=strdup(val);
        else if(!strcasecmp("-simrep",key)) p->simrep=atoi(val);
        else if(!strcasecmp("-is2Dsim",key)) p->is2Dsim=atoi(val);
        else if(!strcasecmp("-numsites",key)) p->numsites=atoi(val);
        else if(!strcasecmp("-RD",key)) p->RD=atof(val);
        else if(!strcasecmp("-tdiv",key)) p->tdiv=atof(val);
        else if(!strcasecmp("-t1",key)) p->t1=atof(val);
        else if(!strcasecmp("-t2",key)) p->t2=atof(val);
        else if(!strcasecmp("-p_inv",key)) p->p_inv=atof(val);
        else if(!strcasecmp("-isthreading",key)) p->isthreading=atoi(val);
        else if(!strcasecmp("-isuchar",key)) p->isuchar=atoi(val);
        else if(!strcasecmp("-dobinary",key)) p->dobinary=atoi(val);
        else if(!strcasecmp("-par",key)){
            stringstream ss(strdup(val));
            while (ss.good() && k<=8){
                string substr;
                getline(ss, substr, ',');
                p->par[k] = stod(substr);
                if (k >= 5){
                    sumpi = sumpi+p->par[k];
                }
                k = k+1;
            }
        }
        else{
            fprintf(stderr,"\t Unknown parameter key:%s val:%s\n",key,val);
            free(p);
            return NULL;
        }
        
        ++argv;
    }
    if ((p->simrep != NULL && p->vcfname != NULL) || (p->simrep == NULL && p->vcfname == NULL)){
        fprintf(stderr,"\t Please specify whether you want to conduct simulation or to analyse vcf file\n");
        free(p);
        return NULL;
    }else if (!strcasecmp(p->model,"GTR") && k!=9){
        fprintf(stderr,"\t In GTR model, 9 parameters should be supplied to the substitution rate matrix \n");
        free(p);
        return NULL;
    }else if (p->vcfname !=NULL && strcasecmp(p->method,"geno")){
        fprintf(stderr,"\t Currenly only distAngsd-geno is implemented for real data analyses!\n");
        free(p);
        return NULL;
    }else if(!strcasecmp(p->model,"GTR")){
        for (int i = 5; i < 9; i++){
            p->par[i] = p->par[i]/sumpi;
        }
    }
    return p;
}
