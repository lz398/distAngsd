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
    //filenames
    p->outname = strdup("distAngsdlog");
    // method can be geno, nuc, RandomSEQ, ConsensusSEQ, AmbiguityGT for simulation
    // method can only be geno for vcf
    p->method = strdup("geno");
    // model can be either JC or GTR
    p->model = strdup("JC");
    p->glfname = NULL;
    p->isthreading = 0;
    
    // Specify the model
    p->is2Dinfer = 0;
    for (int i=0;i<5;i++){
        p->par[i] = 1.0;
    }
    for (int i=5;i<9;i++){
        p->par[i] = 0.25;
    }
    
    // For simulation
    p->simrep = -1;
    p->errorrate = 0.002;
    p->numsites = 1000000;
    p->RD = 1;
    p->tdiv = 1;
    p->t1 = 0.4;
    p->t2 = 0.25;
    p->p_inv = 0;
    
    // For real data
    p->vcfname = NULL;
    p->mpileupname = NULL;
    
    // For table data
    p->tabname = NULL;
    
    
    p->isuchar = 0;
    p->dobinary = 1;
    p->tabuchar = 0;
    p->tabbinary = 1;
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
        else if(!strcasecmp("-outglf",key)) p->glfname=strdup(val);
        else if(!strcasecmp("-vcf",key)) p->vcfname=strdup(val);
        else if(!strcasecmp("-mpileup",key)) p->mpileupname=strdup(val);
        else if(!strcasecmp("-simrep",key)) p->simrep=atoi(val);
        else if(!strcasecmp("-is2Dinfer",key)) p->is2Dinfer=atoi(val);
        else if(!strcasecmp("-e",key)) p->errorrate=atof(val);
        else if(!strcasecmp("-numsites",key)) p->numsites=atoi(val);
        else if(!strcasecmp("-RD",key)) p->RD=atof(val);
        else if(!strcasecmp("-tdiv",key)) p->tdiv=atof(val);
        else if(!strcasecmp("-t1",key)) p->t1=atof(val);
        else if(!strcasecmp("-t2",key)) p->t2=atof(val);
        else if(!strcasecmp("-p_inv",key)) p->p_inv=atof(val);
        else if(!strcasecmp("-isthreading",key)) p->isthreading=atoi(val);
        else if(!strcasecmp("-outuchar",key)) p->isuchar=atoi(val);
        else if(!strcasecmp("-outbin",key)) p->dobinary=atoi(val);
        else if(!strcasecmp("-inglf",key)) p->tabname=strdup(val);
        else if(!strcasecmp("-inuchar",key)) p->tabuchar=atoi(val);
        else if(!strcasecmp("-inbin",key)) p->tabbinary=atoi(val);
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
    int sum = (int)(p->simrep > 0) + (int)(p->vcfname != NULL) + (int)(p->mpileupname != NULL) + (int)(p->tabname !=NULL);
    if (sum !=1){
        fprintf(stderr,"\t Please specify whether you want to conduct simulation or to analyse vcf file.\n");
        free(p);
        return NULL;
    }else if((p->simrep > 0) && (p->tdiv - p->t1 - p->t2 < 0)){
        fprintf(stderr,"\t tdiv should be larger or equal to t1+t2 in simulations.\n");
        free(p);
        return NULL;
    }else if (!strcasecmp(p->model,"GTR") && k!=9){
        fprintf(stderr,"\t In GTR model, 9 parameters should be supplied to the substitution rate matrix.\n");
        free(p);
        return NULL;
    }else if (p->vcfname !=NULL && strcasecmp(p->method,"geno")){
        fprintf(stderr,"\t Currenly only distAngsd-geno is implemented for vcf data analyses!\n");
        free(p);
        return NULL;
    }else if (p->mpileupname !=NULL && strcasecmp(p->method,"nuc")){
        fprintf(stderr,"\t Currenly only distAngsd-nuc is implemented for mpileup data analyses!\n");
        free(p);
        return NULL;
    }else if (p->tabname !=NULL && strcasecmp(p->method,"geno")){
        fprintf(stderr,"\t Currenly only distAngsd-geno is implemented for table glf data analyses!\n");
        free(p);
        return NULL;
    }else if ((p->is2Dinfer==1) && ((strcasecmp(p->model,"GTR")) || (!strcasecmp(p->method,"AmbiguityGT")) || (!strcasecmp(p->method,"NoAmbiguityGT")))){
        fprintf(stderr,"\t Symmetric models, such as JC, are proved to have multi-solutions in 2D inferences, please try GTR.\n \t AmbiguityGT or NoAmbiguity is not suitable for 2D inferences, please try other methods.\n");
        free(p);
        return NULL;
    }else if(!strcasecmp(p->model,"GTR")){
        for (int i = 5; i < 9; i++){
            p->par[i] = p->par[i]/sumpi;
        }
    }
    return p;
}
