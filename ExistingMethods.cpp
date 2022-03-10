#include <iostream>
#include "ExistingMethods.h"
#include "brent.h"
// Ambiguity Inference
//Inference: total likelihood based on the inferred joint genotype distribution + _Ambiguity
double likeGLwithtwoDSFS_Ambiguity(double twoDSFS[10][10], double t, double par[8])
{
    int g1[2], g2[2];
    double like, totlike=0.0;
    //std::cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            if (g1[0] == g1[1] && g2[0] == g2[1]){
                like = log(pi(g1[0])*PMAT(g1[0],g2[0]));
            }else if(g1[0] == g1[1]){
                like = log(pi(g1[0])*(PMAT(g1[0],g2[0])+PMAT(g1[0],g2[1])));
            }else if(g2[0] == g2[1]){
                like = log(pi(g2[0])*(PMAT(g2[0],g1[0])+PMAT(g2[0],g1[1])));
            }else{
                like = log(pi(g1[0])*(PMAT(g1[0],g2[0])+PMAT(g1[0],g2[1]))+pi(g1[1])*(PMAT(g1[1],g2[0])+PMAT(g1[1],g2[1])));
            }
            totlike += like*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}


//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT_Ambiguity(double t)
{
    return likeGLwithtwoDSFS_Ambiguity(GLOBtwoDSFS, t,  GLOBpar);
}

//Inference: Estimation of divergence time t based on joint genotype distribution + _Ambiguity
double estimateT_Ambiguity(double twoDSFS[10][10], double *t, double parameters[])
{
    /*OK this is stupid*/
    for(int i=0;i<10;i++){
        for (int j=0; j<10; j++){
            GLOBtwoDSFS[i][j]=twoDSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT_Ambiguity, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

/* Ambiguity Genotype*/
//Inference: Main EM algorithm for joint genotype distribution
void estimate2DSFS_AmbiguityGT(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites)
{
    double ESFS2[10][10], ptemp, d;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            twoDSFS[i][j]=0.00;
        }
    }
    int sum = 0;
    for (int s=0; s<numsites; s++){
        if (SDATA[s]>0){
            int con_gt1 = 0;
            int con_gt2 = 0;
            int compare_t1 = 0;
            int compare_t2 = 0;
            for (int i=1;i<10;i++){
                if (GLDATA[s][con_gt1]<GLDATA[s][i]){
                    con_gt1 = i;
                    compare_t1 = 0;
                }else if(GLDATA[s][con_gt1]==GLDATA[s][i]){
                    double criti_frac1 = (double)1-(double)1/(double)(compare_t1+2);
                    double u1 = uniform();
                    if (u1 > criti_frac1){
                        con_gt1 = i;
                    }
                    compare_t1 = compare_t1 + 1;
                }
                if (GLDATA[s][10+con_gt2]<GLDATA[s][10+i]){
                    con_gt2 = i;
                    compare_t2 = 0;
                }else if(GLDATA[s][10+con_gt2]==GLDATA[s][10+i]){
                    double criti_frac2 = (double)1-(double)1/(double)(compare_t2+2);
                    double u2 = uniform();
                    if (u2 > criti_frac2){
                        con_gt2 = i;
                    }
                    compare_t2 = compare_t2 + 1;
                }
            }
            //            std::cout<<"individual 1:\n";
            //            for (int i=1;i<10;i++){
            //                std::cout<<GLDATA[s][i]<<"\t";
            //            }
            //
            //            std::cout<<"individual 2:\n";
            //            for (int i=1;i<10;i++){
            //                std::cout<<GLDATA[s][10+i]<<"\t";
            //            }
            //            std::cout<<"\n";
            //            std::cout << s <<" "<< con_gt1 <<" "<<con_gt2<<"\n";
            twoDSFS[con_gt1][con_gt2] = twoDSFS[con_gt1][con_gt2] + 1;
            sum = sum + 1;
        }
    }
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            twoDSFS[i][j]=twoDSFS[i][j]/(double)sum;
        }
    }
}


/* No Ambiguity Genotype*/
//Inference: Main EM algorithm for joint genotype distribution
void estimate2DSFS_NoAmbiguityGT(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites)
{
    double ESFS2[10][10], ptemp, d;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            twoDSFS[i][j]=0.00;
        }
    }
    int sum = 0;
    for (int s=0; s<numsites; s++){
        if (SDATA[s]>0){
            int con_gt1 = 0;
            int con_gt2 = 0;
            int compare_t1 = 0;
            int compare_t2 = 0;
            for (int i=1;i<10;i++){
                if (GLDATA[s][con_gt1]<GLDATA[s][i]){
                    con_gt1 = i;
                    compare_t1 = 0;
                }else if(GLDATA[s][con_gt1]==GLDATA[s][i]){
                    double criti_frac1 = (double)1-(double)1/(double)(compare_t1+2);
                    double u1 = uniform();
                    if (u1 > criti_frac1){
                        con_gt1 = i;
                    }
                    compare_t1 = compare_t1 + 1;
                }
                if (GLDATA[s][10+con_gt2]<GLDATA[s][10+i]){
                    con_gt2 = i;
                    compare_t2 = 0;
                }else if(GLDATA[s][10+con_gt2]==GLDATA[s][10+i]){
                    double criti_frac2 = (double)1-(double)1/(double)(compare_t2+2);
                    double u2 = uniform();
                    if (u2 > criti_frac2){
                        con_gt2 = i;
                    }
                    compare_t2 = compare_t2 + 1;
                }
            }
            int g1[2], g2[2];
            findgenotypes_from_index(con_gt1, g1);
            findgenotypes_from_index(con_gt2, g2);
            if ((g1[0] == g1[1]) && (g2[0] == g2[1])){
                twoDSFS[con_gt1][con_gt2] = twoDSFS[con_gt1][con_gt2] + 1;
                sum = sum + 1;
            }
            
        }
    }
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            twoDSFS[i][j]=twoDSFS[i][j]/(double)sum;
        }
    }
}



/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + _Ambiguity*/
double testtwoDSFS_Ambiguity(double RD, int numsites, double tdiv, double t1, double t2, double errorrate)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10], twoDSFS1[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    pi(0)=0.2184;
    pi(1)=0.2606;
    pi(2)=0.3265;
    pi(3)=0.1946;
    parameters[0]=2.0431;
    parameters[1]=0.0821;
    parameters[2]=0.0000;
    parameters[3]=0.0670;
    parameters[4]=0.0000;
    for (int i=5;i<8; i++)
        parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //Estimate T
    estimateT_Ambiguity(twoDSFS, &t, parameters);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    return t;
}

// Ambiguity Inference
//Inference: total likelihood based on the inferred joint genotype distribution + _Ambiguity
double likeGLwithtwoDSFS_Ambiguity_v2(double twoDSFS[10][10], double t, double par[8])
{
    int g1[2], g2[2];
    double like, totlike=0.0;
    //std::cout<<t<<"\n";
    
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    
    gettransitionprobmatGTR(t);  /*then pi[] has to be initialized somehow else*/
    
    for (int i=0; i<10; i++){
        findgenotypes_from_index(i, g1);
        for (int j=0; j<10; j++){
            findgenotypes_from_index(j, g2);
            if (g1[0] == g1[1] && g2[0] == g2[1]){
                like = log(4.0*pi(g1[0])*PMAT(g1[0],g2[0]));
            }else if(g1[0] == g1[1]){
                like = log(2.0*pi(g1[0])*(PMAT(g1[0],g2[0])+PMAT(g1[0],g2[1])));
            }else if(g2[0] == g2[1]){
                like = log(2.0*pi(g2[0])*(PMAT(g2[0],g1[0])+PMAT(g2[0],g1[1])));
            }else{
                like = log(pi(g1[0])*(PMAT(g1[0],g2[0])+PMAT(g1[0],g2[1]))+pi(g1[1])*(PMAT(g1[1],g2[0])+PMAT(g1[1],g2[1])));
            }
            totlike += like*twoDSFS[i][j];
        }
    }
    //    printf("%lf: %lf\n",t,totlike);
    return -totlike;
}

//Inference: Calculate the likelihood for divergence t, joint genotype distribution
double likelihoodforT_Ambiguity_v2(double t)
{
    return likeGLwithtwoDSFS_Ambiguity_v2(GLOBtwoDSFS, t,  GLOBpar);
}

//Inference: Estimation of divergence time t based on joint genotype distribution + _Ambiguity
double estimateT_Ambiguity_v2(double twoDSFS[10][10], double *t, double parameters[])
{
    /*OK this is stupid*/
    for(int i=0;i<10;i++){
        for (int j=0; j<10; j++){
            GLOBtwoDSFS[i][j]=twoDSFS[i][j];
        }
    }
    for (int i=0; i<8; i++){
        GLOBpar[i]=parameters[i];
    }
    
    double MLV = brent(0.0000001, 0.1, 10.0, likelihoodforT_Ambiguity_v2, 0.000001, t);
    //printf("Max. like value = %lf, t=%lf\n",MLV,*t);
    return MLV;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution + _Ambiguity*/
double testtwoDSFS_Ambiguity_v2(double RD, int numsites, double tdiv, double t1, double t2, double errorrate)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    pi(0)=0.2184;
    pi(1)=0.2606;
    pi(2)=0.3265;
    pi(3)=0.1946;
    parameters[0]=2.0431;
    parameters[1]=0.0821;
    parameters[2]=0.0000;
    parameters[3]=0.0670;
    parameters[4]=0.0000;
    for (int i=5;i<8; i++)
        parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //Estimate T
    estimateT_Ambiguity_v2(twoDSFS, &t, parameters);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    return t;
}

//ngsDist
//Basically, we calculate Fst and use that to estimate the ratio of coalescence time.
void gen_dist(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites, double &dist, double &dist_JC, double &dist1, double &dist1_JC){
    int g1[2], g2[2];
    double score[10][10];
    dist = 0;
    for (int i = 0; i < 10; i++){
        findgenotypes_from_index(i, g1);
        for (int j = 0; j < 10; j++){
            findgenotypes_from_index(j, g2);
            if (g1[0] == g2[0] && g1[1] == g2[1]){
                score[i][j] = 0.0;
            }else if ((g1[0] == g2[0] || g1[1] == g2[1] || g1[0] == g2[1] || g1[1] == g2[0]) && (g1[0] == g1[1] || g2[0] == g2[1])){
                score[i][j] = 0.5;
            }else if (g1[0] == g2[0] || g1[1] == g2[1] || g1[0] == g2[1] || g1[1] == g2[0]){
                score[i][j] = (double) 1.0/ (double) 3.0;
            }else if ((g1[0] == g1[1]) && (g2[0] == g2[1])){
                score[i][j] = 1.0;
            }else if (g1[0] == g1[1] || g2[0] == g2[1]){
                score[i][j] = (double) 3.0/ (double) 4.0;
            }else{
                score[i][j] = 0.5;
            }
            dist += score[i][j]*twoDSFS[i][j];
        }
    }
    dist_JC = -log(1 - (dist * 4/3)) * 3/4;
    dist = -log(1-dist);
    dist1 = 0;
    double dist2;
    double psum;
    for (int s = 0; s<numsites; s++){
        if (SDATA[s]>0){
            dist2 = 0;
            psum = 0;
            for (int i = 0; i < 10; i++){
                for (int j = 0; j < 10; j++){
                    dist2 += score[i][j]*GLDATA[s][i]*GLDATA[s][10+j];
                    psum += GLDATA[s][i]*GLDATA[s][10+j];
                }
            }
            dist1 += dist2/psum;
        }
    }
    dist1 = dist1/(double)eff_numsites;
    dist1_JC = -log(1 - (dist1 * 4/3)) * 3/4;
    dist1 = -log(1-dist1);
}


double testtwoDSFS_ngsDist(double RD, int numsites, double tdiv, double t1, double t2, double errorrate)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    double dist, dist1, dist_JC, dist1_JC;
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    pi(0)=0.2184;
    pi(1)=0.2606;
    pi(2)=0.3265;
    pi(3)=0.1946;
    parameters[0]=2.0431;
    parameters[1]=0.0821;
    parameters[2]=0.0000;
    parameters[3]=0.0670;
    parameters[4]=0.0000;
    for (int i=5;i<8; i++)
        parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //    for (int i=0; i<10; i++){
    //        for (int j=0; j<10; j++){
    //            std::cout << twoDSFS[i][j] <<"\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //    //Estimate T
    //    estimateT_Ambiguity_v2(twoDSFS, &t, parameters);
    //    std::cout<<"Estimated t = "<<t<<"\n";
    //    gen_dist(twoDSFS);
    gen_dist(twoDSFS, GLDATA, SDATA, eff_numsites, numsites, dist, dist_JC, dist1, dist1_JC);
    std::cout << "2DSFS ngsDist result is "<<dist<<"\n";
    double t_dist = (t1+t2)/2/(1-dist);
    std::cout << "2DSFS ngsDist inferred divergence time t is "<<t_dist<<"\n";
    std::cout << "2DSFS ngsDist JC result is "<<dist_JC<<"\n";
    std::cout << "NgsDist result is "<<dist1<<"\n";
    double t_dist1 = (t1+t2)/2/(1-dist1);
    std::cout << "NgsDist inferred divergence time t is "<<t_dist1<<"\n";
    std::cout << "NgsDist JC result is "<<dist1_JC<<"\n";
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    return t;
}

// Another Extension of ngsDist: Basically only account for "SNP"
void gen_dist_snp(double twoDSFS[10][10], double **GLDATA, int *SDATA, int eff_numsites, int numsites, double &dist, double &dist_JC, double &dist1, double &dist1_JC){
    int g1[2], g2[2];
    double score[10][10];
    dist = 0;
    double pcsum = 0;
    for (int i = 0; i < 10; i++){
        findgenotypes_from_index(i, g1);
        for (int j = 0; j < 10; j++){
            findgenotypes_from_index(j, g2);
            if (g1[0] == g2[0] && g1[1] == g2[1] && g1[0] == g1[1]){
                score[i][j] = 0.0;
            }else{
                if(g1[0] == g2[0] && g1[1] == g2[1]){
                    score[i][j] = 0.0;
                }else if ((g1[0] == g2[0] || g1[1] == g2[1] || g1[0] == g2[1] || g1[1] == g2[0]) && (g1[0] == g1[1] || g2[0] == g2[1])){
                    score[i][j] = 0.5;
                }else if (g1[0] == g2[0] || g1[1] == g2[1] || g1[0] == g2[1] || g1[1] == g2[0]){
                    score[i][j] = (double) 1.0/ (double) 3.0;
                }else if ((g1[0] == g1[1]) && (g2[0] == g2[1])){
                    score[i][j] = 1.0;
                }else if (g1[0] == g1[1] || g2[0] == g2[1]){
                    score[i][j] = (double) 3.0/ (double) 4.0;
                }else{
                    score[i][j] = 0.5;
                }
                pcsum = pcsum + twoDSFS[i][j];
            }
            dist += score[i][j]*twoDSFS[i][j];
        }
    }
    //    for (int i=0;i<9;i++){
    //        for (int j=0;j<9;j++){
    //            std::cout<<score[i][j]<<"\t";
    //        }
    //        std::cout<<"\n";
    //    }
    //    std::cout<<pcsum<<"hahaha\n";
    //    std::cout<<dist<<"hahaha\n";
    dist = dist/pcsum;
    dist_JC = -log(1 - (dist * 4/3)) * 3/4;
    //    dist = -log(1-dist);
    dist1 = 0;
    double eff_num = 0;
    double dist2;
    double psum;
    for (int s = 0; s<numsites; s++){
        if (SDATA[s]>0){
            dist2 = 0;
            psum = 0;
            pcsum = 0;
            for (int i = 0; i < 10; i++){
                findgenotypes_from_index(i, g1);
                for (int j = 0; j < 10; j++){
                    findgenotypes_from_index(j, g2);
                    if (g1[0] == g2[0] && g1[1] == g2[1] && g1[0] == g1[1]){
                        psum += GLDATA[s][i]*GLDATA[s][10+j];
                    }else{
                        psum += GLDATA[s][i]*GLDATA[s][10+j];
                        pcsum += GLDATA[s][i]*GLDATA[s][10+j];
                        dist2 += score[i][j]*GLDATA[s][i]*GLDATA[s][10+j];
                    }
                }
            }
            eff_num += pcsum/psum;
            dist1 += dist2/psum;
        }
    }
    dist1 = dist1/eff_num;
    dist1_JC = -log(1 - (dist1 * 4/3)) * 3/4;
    //    dist1 = -log(1-dist1);
}

double testtwoDSFS_ngsDist_snp(double RD, int numsites, double tdiv, double t1, double t2, double errorrate)
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    double dist, dist1, dist_JC, dist1_JC;
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    pi(0)=0.2184;
    pi(1)=0.2606;
    pi(2)=0.3265;
    pi(3)=0.1946;
    parameters[0]=2.0431;
    parameters[1]=0.0821;
    parameters[2]=0.0000;
    parameters[3]=0.0670;
    parameters[4]=0.0000;
    for (int i=5;i<8; i++)
        parameters[i]=pi(i-5);
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_EM(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //    for (int i=0; i<10; i++){
    //        for (int j=0; j<10; j++){
    //            std::cout << twoDSFS[i][j] <<"\t";
    //        }
    //        std::cout << "\n";
    //    }
    
    //    //Estimate T
    //    estimateT_Ambiguity_v2(twoDSFS, &t, parameters);
    //    std::cout<<"Estimated t = "<<t<<"\n";
    //    gen_dist(twoDSFS);
    gen_dist_snp(twoDSFS, GLDATA, SDATA, eff_numsites, numsites, dist, dist_JC, dist1, dist1_JC);
    std::cout << "2DSFS ngsDist result is "<<dist<<"\n";
    double t_dist = (t1+t2)/2/(1-dist);
    std::cout << "2DSFS ngsDist inferred divergence time t is "<<t_dist<<"\n";
    std::cout << "2DSFS ngsDist JC result is "<<dist_JC<<"\n";
    std::cout << "NgsDist result is "<<dist1<<"\n";
    double t_dist1 = (t1+t2)/2/(1-dist1);
    std::cout << "NgsDist inferred divergence time t is "<<t_dist1<<"\n";
    std::cout << "NgsDist JC result is "<<dist1_JC<<"\n";
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    return t;
}


/*Consensus Read Inference*/
void simSEQs_reads_consensus(double **SEQDATA, int genotypes[2], int site, int species, double e, double RD)
{
    int ReadDepth = ReaddepthGenerator(RD, 0);
    int nuccount[4]={0,0,0,0};
    int out;
    for(int k=0;k<4;k++){
        SEQDATA[site][4*species+k] = 0;
    }
    for(int i=0;i<ReadDepth;i++){
        out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e,0);
        //        double e1 = -e*log(1-uniform());
        //        if (e1 > 0.5) {e1=0.5;}
        //        if (i!=j){
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,0);
        //        }else{
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,1);
        //        }
        nuccount[out] = nuccount[out]+1;
    }
    if (ReadDepth>0){
        int con_out = 0;
        int compare_t = 0;
        for (int i=1;i<4;i++){
            if (nuccount[con_out]<nuccount[i]){
                con_out = i;
                compare_t = 0;
            }else if(nuccount[con_out]==nuccount[i]){
                double criti_frac = (double)1-(double)1/(double)(compare_t+2);
                double u = uniform();
                if (u > criti_frac){
                    con_out = i;
                }
                compare_t = compare_t + 1;
            }
        }
        SEQDATA[site][4*species+con_out] = 1;
    }
    //    double a = log(1.0-e);
    //    double b = log(e)-log(3.0);
    //    for(int k=0;k<4;k++){
    //        if (k == con_out){
    //            SEQDATA[site][4*species+k] = SEQDATA[site][4*species+k] + a;
    //        }else{
    //            SEQDATA[site][4*species+k] = SEQDATA[site][4*species+k] + b;
    //        }
    //    }
    //
    //    if (ReadDepth>0){
    //        for(int k=0;k<4;k++){
    //            SEQDATA[site][4*species+k] = exp(SEQDATA[site][4*species+k]);
    //            //std::cout<<SEQDATA[site][4*species+k]<<"\t";
    //        }
    //    }
    //std::cout<<"\n";
}

void simSEQs_reads_random(double **SEQDATA, int genotypes[2], int site, int species, double e, double RD)
{
    int ReadDepth = ReaddepthGenerator(RD, 0);
    int nuccount[4]={0,0,0,0};
    int out;
    for(int k=0;k<4;k++){
        SEQDATA[site][4*species+k] = 0;
    }
    for(int i=0;i<ReadDepth;i++){
        out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e,0);
        //        double e1 = -e*log(1-uniform());
        //        if (e1 > 0.5) {e1=0.5;}
        //        if (i!=j){
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,0);
        //        }else{
        //            out = simSEQs_per_read_v1(SEQDATA, genotypes, site, species, e1,1);
        //        }
        nuccount[out] = nuccount[out]+1;
    }
    if (ReadDepth>0){
        double p_cumnuc=(double)nuccount[0]/(double)ReadDepth;
        int rand_out = 0;
        double u = uniform();
        while (u > p_cumnuc){
            rand_out = rand_out + 1;
            p_cumnuc = p_cumnuc + (double)nuccount[rand_out]/(double)ReadDepth;
        }
        SEQDATA[site][4*species+rand_out] = 1;
    }
}


void simulateGLsTwoSpeciesSEQ_random(double RD, int numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (int i=0; i<numsites; i++)
        ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (int i=0; i<numsites; i++){
        simnucleotides(ancDATA[i], SIMMAT);
        //printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    std::cout<<"Ancient Data was derived!\n";
    
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    // makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for (int i=0; i<numsites; i++){
        simpoly(ancDATA[i][0], SIMMAT, genotypes);
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        // simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_random(SEQDATA, genotypes, i, 0, errorrate, RD);
        
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (int i=0; i<numsites; i++){
        simpoly(ancDATA[i][1], SIMMAT, genotypes);
        //simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_random(SEQDATA, genotypes, i, 1, errorrate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    /*then we free memory*/
    for (int i=0; i<numsites; i++)
        free(ancDATA[i]);
    free(ancDATA);
}

void simulateGLsTwoSpeciesSEQ_consensus(double RD, int numsites, double errorrate, double tdiv, double t1, double t2, double **SEQDATA, double (*pijt)(double t, double *par, int from, int to), double *par)
{
    double  SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int **ancDATA, genotypes[2];
    
    if (tdiv<t1/2.0 || tdiv<t2/2.0) {printf("Divergence time misspecified"); exit(-1);}
    
    /*first we allocate memory*/
    ancDATA = (int **) malloc(numsites*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (int i=0; i<numsites; i++)
        ancDATA[i]=(int *) malloc(2*(sizeof(int)));
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(tdiv-t1-t2); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (int i=0; i<numsites; i++){
        simnucleotides(ancDATA[i], SIMMAT);
        //printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    std::cout<<"Ancient Data was derived!\n";
    
    /*k=0;
     for (i=0; i<numsites; i++){
     if (ancDATA[i][0] != ancDATA[i][1]) k++;
     }
     printf("ancestral numdif/length = %lf (%i)\n",(double)k/(double)numsites, k);?*/
    
    
    /*then we simulate the genotype in species1 and the resulting GL data*/
    gettransitionprobmatGTR(t1);
    makepolyMAT(SIMMAT);
    // makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for (int i=0; i<numsites; i++){
        simpoly(ancDATA[i][0], SIMMAT, genotypes);
        //    if (ancDATA[i][0] != genotypes[0]) k++;
        //    if (ancDATA[i][0] != genotypes[1]) k++;
        // simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 0, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_consensus(SEQDATA, genotypes, i, 0, errorrate, RD);
        
        //printf(" site %i, genotypes species 1: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    
    //printf("%lf differences were added on branch 1\n",(double)k/(2.0*(double)numsites));
    
    
    /*then we simulate the genotype in species2 and the resulting GL data*/
    gettransitionprobmatGTR(t2);
    makepolyMAT(SIMMAT);
    for (int i=0; i<numsites; i++){
        simpoly(ancDATA[i][1], SIMMAT, genotypes);
        //simSEQs(SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), i, 1, switchmatrix, simswitchmatrix, errorrate);
        simSEQs_reads_consensus(SEQDATA, genotypes, i, 1, errorrate, RD);
        //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
    }
    /*then we free memory*/
    for (int i=0; i<numsites; i++)
        free(ancDATA[i]);
    free(ancDATA);
}

double testsimSEQDATA_random(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[8])
{
    double t, MLV, parameters[8];
    
    
    //SetSeed(6);
    SetSeed(rand()%30000+1);
    //    for (int i=0; i<4; i++)
    //    pi(i)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    //    pi(0)=0.2184;
    //    pi(1)=0.2606;
    //    pi(2)=0.3265;
    //    pi(3)=0.1946;
    //    parameters[0]=2.0431;
    //    parameters[1]=0.0821;
    //    parameters[2]=0.0000;
    //    parameters[3]=0.0670;
    //    parameters[4]=0.0000;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    
    
    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    SEQDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        SEQDATA[i] = (double *) malloc(8 * sizeof(double));
    
    SEQ_SDATA = (int *) malloc(numsites * sizeof(int));
    
    simulateGLsTwoSpeciesSEQ_random(RD, numsites, errorrate,  tdiv,  t1,  t2, SEQDATA, pijtGTR, parameters);
    
    int eff_numsites = FilterSitesSEQ(SEQDATA, SEQ_SDATA, numsites);
    
    for (int i=0; i<8; i++)
        GLOBpar[i]=parameters[i];
    globnumsites=numsites;
    globerror=errorrate;
    MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ_v1, 0.000001, &t);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    for (int i = 0; i < numsites; i++)
        free(SEQDATA[i]);
    free(SEQDATA);
    free(SEQ_SDATA);
    return t;
}

double testsimSEQDATA_consensus(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9])
{
    double t, MLV, parameters[8];
    
    
    //SetSeed(6);
    SetSeed(rand()%30000+1);
    //    for (int i=0; i<4; i++)
    //    pi(i)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    //    pi(0)=0.2184;
    //    pi(1)=0.2606;
    //    pi(2)=0.3265;
    //    pi(3)=0.1946;
    //    parameters[0]=2.0431;
    //    parameters[1]=0.0821;
    //    parameters[2]=0.0000;
    //    parameters[3]=0.0670;
    //    parameters[4]=0.0000;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    
    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    SEQDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        SEQDATA[i] = (double *) malloc(8 * sizeof(double));
    
    SEQ_SDATA = (int *) malloc(numsites * sizeof(int));
    
    simulateGLsTwoSpeciesSEQ_consensus(RD, numsites, errorrate,  tdiv,  t1,  t2, SEQDATA, pijtGTR, parameters);
    
    int eff_numsites = FilterSitesSEQ(SEQDATA, SEQ_SDATA, numsites);
    
    for (int i=0; i<8; i++)
        GLOBpar[i]=parameters[i];
    globnumsites=numsites;
    globerror=errorrate;
    MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ_v1, 0.000001, &t);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    for (int i = 0; i < numsites; i++)
        free(SEQDATA[i]);
    free(SEQDATA);
    free(SEQ_SDATA);
    return t;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution*/
double testtwoDSFS_AmbiguityGT(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9])
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    //    pi(0)=0.2184;
    //    pi(1)=0.2606;
    //    pi(2)=0.3265;
    //    pi(3)=0.1946;
    //    parameters[0]=2.0431;
    //    parameters[1]=0.0821;
    //    parameters[2]=0.0000;
    //    parameters[3]=0.0670;
    //    parameters[4]=0.0000;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    //
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_AmbiguityGT(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //Estimate T
    estimateT_Ambiguity(twoDSFS, &t, parameters);
    //estimateT(twoDSFS, &t, parameters);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    
    free(SDATA);
    return t;
}

/*Simulation + Inference: Simulation and estimation of divergence time t based on joint genotype distribution*/
double testtwoDSFS_NoAmbiguityGT(double RD, int numsites, double tdiv, double t1, double t2, double errorrate, double par[9])
{
    double **GLDATA, t, parameters[8], twoDSFS[10][10];
    int *SDATA;
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
        GLDATA[i] =(double *) malloc(20 * sizeof(double));
    
    SDATA = (int *) malloc(numsites * sizeof(int));
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    //    pi(0)=0.25;
    //    pi(1)=0.25;
    //    pi(2)=0.25;
    //    pi(3)=0.25;
    //    for (int i=0; i<5; i++)
    //    parameters[i]=1.0 /*+ (double)i*/;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    //    pi(0)=0.2184;
    //    pi(1)=0.2606;
    //    pi(2)=0.3265;
    //    pi(3)=0.1946;
    //    parameters[0]=2.0431;
    //    parameters[1]=0.0821;
    //    parameters[2]=0.0000;
    //    parameters[3]=0.0670;
    //    parameters[4]=0.0000;
    //    for (int i=5;i<8; i++)
    //    parameters[i]=pi(i-5);
    for (int i=5;i<9;i++){
        pi[i-5] = par[i];
    }
    for (int i=0;i<8;i++){
        parameters[i] = par[i];
    }
    
    //simulate data
    simulateGLsTwoSpecies(RD, numsites, errorrate,  tdiv,  t1,  t2, GLDATA, pijtGTR, parameters);
    
    //    for (int i = 0; i < numsites; i++){
    //        for (int j = 0; j < 20; j++){
    //            std::cout << GLDATA[i][j] << "\t";
    //        }
    //        std::cout << "\n";
    //    }
    //
    
    //Filter the effective numsites
    int eff_numsites = FilterSites(GLDATA, SDATA, numsites);
    
    //Estimate 2DSFS in 10x10 matrix
    estimate2DSFS_NoAmbiguityGT(twoDSFS, GLDATA, SDATA, eff_numsites, numsites);
    
    //Estimate T
    estimateT_Ambiguity(twoDSFS, &t, parameters);
    //estimateT(twoDSFS, &t, parameters);
    std::cout<<"Estimated t = "<<t<<"\n";
    
    
    for (int i = 0; i < numsites; i++)
        free(GLDATA[i]);
    free(GLDATA);
    
    free(SDATA);
    return t;
}
