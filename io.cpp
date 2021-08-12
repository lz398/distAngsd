#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <cmath>
#include <limits>
#include <string>
#include <pthread.h>
#include <cassert>
#include "special.h"

/*we hardcode this to avoid having to do any arithmetic in the conversion and to make sure homozygous types are checked first*/
/*this function is optimized to be fast and not to be pretty*/
/*
 AA=00: 0
 AC=01: 1
 AG=02: 2
 AT=03: 3
 CC=11: 4
 CG=12: 5
 CT=13: 6
 GG=22: 7
 GT=23: 8
 TT=33: 9
 */
int findgenotypeindex(int i, int j)
{
    int k;

    if (i>j){
        k=i;
        i=j;
        j=k;
    }
    if (i==0){
        if (j==0) return 0;
        else if (j==1) return 1;
        else if (j==2) return 2;
        else if (j==3) return 3;
        else {printf("error in genotype index conversion table"); exit(-1);}
    }
    else if (i==1){
        if (j==1) return 4;
        else if (j==2) return 5;
        else if (j==3) return 6;
        else {printf("error in genotype index conversion table"); exit(-1);}
    }
    else if (i==2){
        if (j==2) return 7;
        else if (j==3) return 8;
        else {printf("error in genotype index conversion table"); exit(-1);}
    }
    else if (i==3 && j==3) return 9;
    else {printf("error in genotype index conversion table"); exit(-1);}
}

void findgenotypes_from_index(int inde, int genotype[2])
{
    switch(inde) {
        case 0: genotype[0]=0; genotype[1]=0; break;
        case 1: genotype[0]=0; genotype[1]=1; break;
        case 2: genotype[0]=0; genotype[1]=2; break;
        case 3: genotype[0]=0; genotype[1]=3; break;
        case 4: genotype[0]=1; genotype[1]=1; break;
        case 5: genotype[0]=1; genotype[1]=2; break;
        case 6: genotype[0]=1; genotype[1]=3; break;
        case 7: genotype[0]=2; genotype[1]=2; break;
        case 8: genotype[0]=2; genotype[1]=3; break;
        case 9: genotype[0]=3; genotype[1]=3; break;
        default: {printf("error 2 in genotype index conversion table (case=%i)",inde); exit(-1);}
    }
}
