#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
typedef Eigen::Matrix<double, 10, 10> Matrix10d;
typedef Eigen::Matrix<double, 10, 1> Vector10d;
using namespace std;
using namespace Eigen;
#include "Treetest.h"

const int numsites = 1000000;

string ReadTree(char * filename)
{
    ifstream ifile(filename);
    ostringstream buf;
    char ch;
    char ch1 = ' ';
    char ch2 = ';';
    while(buf&&ifile.get(ch))
        if(ch!=ch1 && ch!=ch2){
            buf.put(ch);
        }
    return buf.str();
}

struct Simleaf{
    string leafname;
    int haplotype[numsites];
    double distance;
    double **GLDATA;
    int *SEQDATA;
    double t;
};
struct SimTree{
    vector<Simleaf> Simleaves;
    vector<string> leafnames;
};

int simhap(int nuc, double SIMMAT[4][4])
{
    double u;
    u=uniform();
    int j=0;
    while (u>SIMMAT[nuc][j]){
        j++;
        if (j>3) {
            printf("Numerical error 4 in simulation algorithm"); exit(-1);
        }
    }
    return j;
}

int treeviewing(char* treefile){
    int m = -1;
    string str = ReadTree(treefile);
    vector<string > nodenames;
    int flg[str.size()];
    flg[0]=0;
    if(str[0]=='('){
        for(int i=1;i<str.size();i++){
            //cout<<str[i]<<"\n";
            if(str[i-1]=='('){
                flg[i]=flg[i-1]+1;
            }else if(str[i]==')'){
                flg[i]=flg[i-1]-1;
            }else{
                flg[i]=flg[i-1];
            }
        }
    }
    int i=0;
    int l=0;
    while(i<str.size()){
        int pt0 = i;
        while((str[i]!=':'&&flg[i]==1)||flg[i]>1){
            i=i+1;
        }
        int pt1 = i;
        if (pt1-pt0>0){
            //cout<<str.substr(pt0,pt1-pt0);
            if(l==0){
                nodenames.push_back(str.substr(pt0,pt1-pt0));
                l=1;
            }
            i=i+1;
        }
        //cout<<"\n";
        pt0=i;
        while(str[i]!=','&&flg[i]==1){
            i=i+1;
        }
        pt1 = i;
        if (pt1-pt0>0){
            //cout<<str.substr(pt0,pt1-pt0);
            if (l==1){
                l=0;
            }
        }
        i=i+1;
    }
    
    int k=0;
    while(k<nodenames.size()){
        if(nodenames[k].size()>1){
            string str=nodenames[k];
            int flg[str.size()];
            flg[0]=0;
            if(str[0]=='('){
                for(int i=1;i<str.size();i++){
                    //cout<<str[i]<<"\n";
                    if(str[i-1]=='('){
                        flg[i]=flg[i-1]+1;
                    }else if(str[i]==')'){
                        flg[i]=flg[i-1]-1;
                    }else{
                        flg[i]=flg[i-1];
                    }
                }
            }
            int i=0;
            int l=0;
            while(i<str.size()){
                int pt0 = i;
                while((str[i]!=':'&&flg[i]==1)||flg[i]>1){
                    i=i+1;
                }
                int pt1 = i;
                if (pt1-pt0>0){
                    if(l==0){
                        nodenames.push_back(str.substr(pt0,pt1-pt0));
                        l=1;
                    }
                    i=i+1;
                }
                pt0=i;
                while(str[i]!=','&&flg[i]==1){
                    i=i+1;
                }
                pt1 = i;
                if (pt1-pt0>0){
                    if (l==1){
                        l=0;
                    }
                }
                i=i+1;
            }
            nodenames.erase(nodenames.begin()+k);
            k=k-1;
        }
        k=k+1;
    }
    return nodenames.size();
}

void treebuilding(string str,SimTree &mytree,double par[]){
    double  SIMMAT[4][4];
    int *splitDATA, hap;
    splitDATA=(int *) malloc(2*(sizeof(int)));
    int hiddenDATA[numsites];
    int flg[str.size()];
    flg[0]=0;
    if(str[0]=='('){
        for(int i=1;i<str.size();i++){
            //cout<<str[i]<<"\n";
            if(str[i-1]=='('){
                flg[i]=flg[i-1]+1;
            }else if(str[i]==')'){
                flg[i]=flg[i-1]-1;
            }else{
                flg[i]=flg[i-1];
            }
        }
    }
    int i=0;
    int l=0;
    Simleaf leaf;
    while(i<str.size()){
        int pt0 = i;
        while((str[i]!=':'&&flg[i]==1)||flg[i]>1){
            i=i+1;
        }
        int pt1 = i;
        if (pt1-pt0>0){
            //cout<<str.substr(pt0,pt1-pt0);
            if(l==0){
                leaf.leafname = str.substr(pt0,pt1-pt0);
                l=1;
            }
            i=i+1;
        }
        //cout<<"\n";
        pt0=i;
        while(str[i]!=','&&flg[i]==1){
            i=i+1;
        }
        pt1 = i;
        if (pt1-pt0>0){
            //cout<<str.substr(pt0,pt1-pt0);
            if (l==1){
                leaf.distance=stod(str.substr(pt0,pt1-pt0));
                mytree.Simleaves.push_back(leaf);
                mytree.leafnames.push_back(leaf.leafname);
                l=0;
            }
        }
        i=i+1;
    }
    cout<<mytree.Simleaves[0].distance<<"\n";
    
    /*then we simulate the two root ancestors for each species*/
    diagonalizeGTR(par); /*note that this code needs to be moved out of the function if JC or HKY functionality is used*/
    gettransitionprobmatGTR(mytree.Simleaves[0].distance); /*instead pi needs to be defined appropriately*/
    makeSIMMAT(SIMMAT);
    for (int i=0; i<numsites; i++){
        simnucleotides(splitDATA, SIMMAT);
        mytree.Simleaves[0].haplotype[i]=splitDATA[0];
        hiddenDATA[i]=splitDATA[1];
        // printf(" site %i, anc. nucs: %i %i\n",i,ancDATA[i][0],ancDATA[i][1]);
    }
    cout<<hiddenDATA[1]<<"HAHAHAH\n";
    for(int j=1;j<3;j++){
        gettransitionprobmatGTR(mytree.Simleaves[j].distance);
        makepolyMAT(SIMMAT);
        for (int i=0; i<numsites; i++){
            mytree.Simleaves[j].haplotype[i]=simhap(hiddenDATA[i], SIMMAT);
        }
    }
    
    for(int j=0;j<mytree.Simleaves.size();j++){
        cout<<"Node Name: "<<mytree.Simleaves[j].leafname<<"\n";
        for (int k=0;k<100;k++){
            cout<<mytree.Simleaves[j].haplotype[k];
        }
        cout<<"\n";
    }
    cout<<"Ancient Data was derived!\n";
}
//void treesprouting(string str,SimTree mytree){
void treesprouting(SimTree &mytree){
    double  SIMMAT[4][4];
    int k=0;
    Simleaf leaf;
    vector<Simleaf > leaves;
    while(k<mytree.Simleaves.size()){
        if(mytree.Simleaves[k].leafname.size()>1){
            cout<<mytree.Simleaves[k].leafname<<"\n";
            //treesprouting(mytree.Simleaves[k].leafname);
            string str=mytree.Simleaves[k].leafname;
            cout<< str<<"\n";
            int flg[str.size()];
            flg[0]=0;
            if(str[0]=='('){
                for(int i=1;i<str.size();i++){
                    //cout<<str[i]<<"\n";
                    if(str[i-1]=='('){
                        flg[i]=flg[i-1]+1;
                    }else if(str[i]==')'){
                        flg[i]=flg[i-1]-1;
                    }else{
                        flg[i]=flg[i-1];
                    }
                }
            }
            int i=0;
            int l=0;
            
            while(i<str.size()){
                int pt0 = i;
                while((str[i]!=':'&&flg[i]==1)||flg[i]>1){
                    i=i+1;
                }
                int pt1 = i;
                if (pt1-pt0>0){
                    if(l==0){
                        leaf.leafname = str.substr(pt0,pt1-pt0);
                        l=1;
                    }
                    i=i+1;
                }
                pt0=i;
                while(str[i]!=','&&flg[i]==1){
                    cout<<str[i];
                    i=i+1;
                }
                pt1 = i;
                cout<<"\n";
                if (pt1-pt0>0){
                    cout<<str.substr(pt0,pt1-pt0);
                    if (l==1){
                        leaf.distance=stod(str.substr(pt0,pt1-pt0));
                        gettransitionprobmatGTR(leaf.distance);
                        makepolyMAT(SIMMAT);
                        for (int i=0; i<numsites; i++){
                            leaf.haplotype[i]=simhap(mytree.Simleaves[k].haplotype[i], SIMMAT);
                        }
                        leaves.push_back(leaf);
                        l=0;
                    }
                }
                i=i+1;
            }
            mytree.Simleaves.erase(mytree.Simleaves.begin()+k);
            mytree.leafnames.erase(mytree.leafnames.begin()+k);
            k=k-1;
        }
        for (int j=0;j<leaves.size();j++){
            mytree.Simleaves.push_back(leaves[j]);
            mytree.leafnames.push_back(leaves[j].leafname);
        }
        leaves.clear();
        k=k+1;
    }
    for (int i=0;i<mytree.Simleaves.size();i++){
        cout<<mytree.Simleaves[i].leafname<<"\n";
        cout<<mytree.Simleaves[i].distance<<"\n";
        for (int j=0;j<100;j++){
            cout<<mytree.Simleaves[i].haplotype[j];
        }
        cout<<"\n";
    }
}

void leavesorting(SimTree &mytree){
    vector<int > leafnames;
    for(int i=0;i<mytree.leafnames.size();i++){
        leafnames.push_back(stoi(mytree.leafnames[i]));
    }
    for(int j=1;j<leafnames.size()-1;j++){
        for(int i=0;i<leafnames.size()-j;i++){
            if(leafnames[i]>leafnames[i+1]){
                iter_swap(leafnames.begin() + i, leafnames.begin() + i+1);
                iter_swap(mytree.leafnames.begin() + i, mytree.leafnames.begin() + i+1);
                iter_swap(mytree.Simleaves.begin() + i, mytree.Simleaves.begin() + i+1);
            }
        }
    }
    //    cout<<"\n";
    //    for(int i=0;i<mytree.leafnames.size();i++){
    //        cout<<leafnames[i]<<"\t"<<mytree.leafnames[i]<<"\t"<<mytree.Simleaves[i].leafname<<"\n";
    //    }
}

void ReadSplitTimes(char * filename,SimTree &mytree){
    ifstream ifs(filename);
    string str;
    int count = 0;
    while (ifs >> str)
    {
        mytree.Simleaves[count].t=stod(str);
        count++;
    }
    ifs.close();
}

void addtreediploid(SimTree& mytree, char* timefile, double errorrate){
    ReadSplitTimes(timefile,mytree);
    double SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int genotypes[2];
    makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for(int i=0;i<mytree.Simleaves.size();i++){
        gettransitionprobmatGTR(mytree.Simleaves[i].t);
        mytree.Simleaves[i].GLDATA = (double **) malloc(numsites * sizeof(double *));
        for (int k = 0; k < numsites; k++)
        mytree.Simleaves[i].GLDATA[k] =(double *) malloc(10 * sizeof(double));
        makepolyMAT(SIMMAT);
        for (int j=0; j<numsites; j++){
            simpoly(mytree.Simleaves[i].haplotype[j], SIMMAT, genotypes);
            simGLs(mytree.Simleaves[i].GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), j, 0, switchmatrix, simswitchmatrix);
            //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
        }
    }
}

void addtreehaploid(SimTree& mytree, char* timefile, double errorrate){
    ReadSplitTimes(timefile,mytree);
    double SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int genotypes[2];
    makeGLswitchmatrix(errorrate, switchmatrix, simswitchmatrix);
    for(int i=0;i<mytree.Simleaves.size();i++){
        gettransitionprobmatGTR(mytree.Simleaves[i].t);
        mytree.Simleaves[i].SEQDATA = (int *) malloc(numsites * sizeof(int));
        makepolyMAT(SIMMAT);
        for (int j=0; j<numsites; j++){
            simpoly(mytree.Simleaves[i].haplotype[j], SIMMAT, genotypes);
            simSEQs(mytree.Simleaves[i].SEQDATA, findgenotypeindex(genotypes[0], genotypes[1]), j, switchmatrix, simswitchmatrix, errorrate);
            //simGLs(mytree.Simleaves[i].GLDATA, findgenotypeindex(genotypes[0], genotypes[1]), j, 0, switchmatrix, simswitchmatrix);
            //printf(" site %i, genotypes species 2: %i %i\n",i,genotypes[0], genotypes[1]);
        }
    }
}

void testsimSEQDATA_tree(double **T,double errorrate, char* treefile, char* timefile){
    double t, parameters[8], SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10], MLV;
    int genotypes[2];
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    SimTree mytree;
    string str = ReadTree(treefile);
    treebuilding(str,mytree,parameters);
    treesprouting(mytree);
    leavesorting(mytree);
    addtreehaploid(mytree,timefile,errorrate);
    
    //This codes tests the program if sampling a single nucleotide
    //printf("T=%.2lf, (t1=%.2lf, t2=%.2lf), e=%.2lf, numb.sites=%i: ",tdiv,t1,t2,errorrate, numsites);
    SEQDATA = (int **) malloc(numsites * sizeof(int *));
    for (int i = 0; i < numsites; i++)
    SEQDATA[i] = (int *) malloc(2 * sizeof(int));
    
    for (int i=0; i<8; i++)
    GLOBpar[i]=parameters[i];
    globnumsites=numsites;
    globerror=errorrate;
    
   
    for(int i=0;i<mytree.Simleaves.size();i++){
        for(int j=0;j<i+1;j++){
            for(int k=0;k<numsites;k++){
                    SEQDATA[k][0] = mytree.Simleaves[i].SEQDATA[k];
                    SEQDATA[k][1] = mytree.Simleaves[j].SEQDATA[k];
            }
            MLV = brent(0.0000001, 0.1, 10.0, likelihoodforTSEQ, 0.000001, &t);
            T[i][j] = t;
            T[j][i] = t;
        }
    }
    mytree.Simleaves.clear();
    mytree.leafnames.clear();
    for (int i = 0; i < numsites; i++)
    free(SEQDATA[i]);
    free(SEQDATA);
}

void testtwoDSFS_tree(double **T,double errorrate, char* treefile, char* timefile){
    double **GLDATA, t, twoDSFS[10][10], parameters[8], SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int genotypes[2];
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    SimTree mytree;
    string str = ReadTree(treefile);
    treebuilding(str,mytree,parameters);
    treesprouting(mytree);
    leavesorting(mytree);
    addtreediploid(mytree,timefile,errorrate);
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));
    for(int i=0;i<mytree.Simleaves.size();i++){
        for(int j=0;j<i+1;j++){
            for(int k=0;k<numsites;k++){
                for(int l=0;l<10;l++){
                    GLDATA[k][l] = mytree.Simleaves[i].GLDATA[k][l];
                    GLDATA[k][10+l] = mytree.Simleaves[j].GLDATA[k][l];
                }
            }
            //Estimate 2DSFS in 10x10 matrix
            estimate2DSFS_EM(twoDSFS, GLDATA, numsites);
            //Estimate T
            estimateT(twoDSFS, &t, parameters);
            T[i][j] = t;
            T[j][i] = t;
        }
    }
    mytree.Simleaves.clear();
    mytree.leafnames.clear();
    for (int i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);
}

void testtwoDSFS_tree_m(double **T,double errorrate, char* treefile, char* timefile){
    double **GLDATA, t, twoDSFS[10][10], parameters[8], SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int genotypes[2];
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    SimTree mytree;
    string str = ReadTree(treefile);
    treebuilding(str,mytree,parameters);
    treesprouting(mytree);
    leavesorting(mytree);
    addtreediploid(mytree,timefile,errorrate);
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));
    for(int i=0;i<mytree.Simleaves.size();i++){
        for(int j=0;j<i+1;j++){
            for(int k=0;k<numsites;k++){
                for(int l=0;l<10;l++){
                    GLDATA[k][l] = mytree.Simleaves[i].GLDATA[k][l];
                    GLDATA[k][10+l] = mytree.Simleaves[j].GLDATA[k][l];
                }
            }
            //Estimate 2DSFS in 10x10 matrix
            estimate2DSFS_EM(twoDSFS, GLDATA, numsites);
            //Estimate T
            estimateT_m(twoDSFS, &t, parameters,mytree.Simleaves[i].t,mytree.Simleaves[j].t);
            T[i][j] = t;
            T[j][i] = t;
        }
    }
    mytree.Simleaves.clear();
    mytree.leafnames.clear();
    for (int i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);
}

void testtwoDSFS_tree_EM(double **T,double errorrate, char* treefile, char* timefile){
    double **GLDATA, t, twoDSFS[10][10], parameters[8], SIMMAT[4][4], switchmatrix[10][10], simswitchmatrix[10][10];
    int genotypes[2];
    pi(0)=0.25;
    pi(1)=0.25;
    pi(2)=0.25;
    pi(3)=0.25;
    for (int i=0; i<5; i++)
    parameters[i]=1.0 /*+ (double)i*/;
    for (int i=5;i<8; i++)
    parameters[i]=pi(i-5);
    
    //SetSeed(666);
    SetSeed(rand()%30000+1);
    
    SimTree mytree;
    string str = ReadTree(treefile);
    treebuilding(str,mytree,parameters);
    treesprouting(mytree);
    leavesorting(mytree);
    addtreediploid(mytree,timefile,errorrate);
    
    //initialize and set simulation parameters
    GLDATA = (double **) malloc(numsites * sizeof(double *));
    for (int i = 0; i < numsites; i++)
    GLDATA[i] =(double *) malloc(20 * sizeof(double));

    GLOBPsum = (double**) malloc(40 * sizeof(double ***));
    for (int i = 0; i < 40; i++){
        GLOBPsum[i] = (double *) malloc(40 * sizeof(double));}
    
    for(int i=0;i<mytree.Simleaves.size();i++){
        for(int j=0;j<i+1;j++){
            for(int k=0;k<numsites;k++){
                for(int l=0;l<10;l++){
                    GLDATA[k][l] = mytree.Simleaves[i].GLDATA[k][l];
                    GLDATA[k][10+l] = mytree.Simleaves[j].GLDATA[k][l];
                }
            }
            double tt[3];
            EMlikelihoodforT(GLDATA, numsites, 2, 0.4, 0.2,&tt[0]);
            T[i][j]=tt[0]+tt[1]+tt[2];
            T[j][i]=T[i][j];
        }
    }
    mytree.Simleaves.clear();
    mytree.leafnames.clear();
    for (int i = 0; i < numsites; i++)
    free(GLDATA[i]);
    free(GLDATA);
    for (int i = 0; i <40; i++)
    free(GLOBPsum[i]);
    free(GLOBPsum);
}
