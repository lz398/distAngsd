#a [![test](https://github.com/lz398/distAngsd/actions/workflows/make.yml/badge.svg)](https://github.com/lz398/distAngsd/actions/workflows/make.yml) distAngsd
# This is a test version, and please let us know if there is any problem while running the code. Thanks a lot!
## Descriptions
<strong>distAngsd</strong> is a software to simulate and infer phylogenetic relationships between two individuals, in which two novel methods (i.e., geno and nuc) are proposed. A bunch of other methods are also implemented, e.g., RandomSEQ, ConsensusSEQ, AmbiguityGT and NoAmbiguity for comparisons. The software can both conduct simulation and analyses real vcf/bcf file given reliable genotype likelihoods are calculated.

The software can do 1-d simulation and inference (only genetic distance <em>t</em> estimation), 2-d simulation and inference (genetic distance <em>t</em> and invariable sites proportion <em>p_inv</em> joint estimation), vcf file read and inference, genotype likelihood table read and inference. Both JC69 are GTR models can be considered. Genotype likelihood and inference logs can be recored either in the format of txt.gz or bin.

The code currently contains Makefile, io.cpp, bfgs.cpp, GLtest.cpp, GL2Dtest.cpp, ExistingMethods.cpp, vcftest.cpp, mpileup.cpp, GL-Reads.cpp, io.h, bfgs.h, GLtest.h, GL2Dtest.h, ExistingMethods.h vcftest.h, mpileup.h, and GL-Reads.h.

## Compilation and Run
### Compilation
make distAngsd

or

make EIGEN=DIRECTORY_TO_YOUR_eigen3 HTSSRC=DIRECTORY_TO_YOUR_htslib.a distAngsd
### Run
./distAngsd -o -method -model -inglf -outglf -vcf -mpileup -simrep -is2Dinfer -isex -p_inv -isthreading -inbin -inuchar -outbin -outuchar -numsites -RD -e -tdiv -t1 -t2 -par

One and only one of the following three options must be provided to specify whether to infer based on simulation or vcf file.
* -vcf specifies the directory of input vcf files.

* -mileup specifies the directory of input mileup bam files.

* -simrep specifies the number of simulation replicates.

* -inglf specifies the directory of input genotype likelihood file. 

* -o specifies the recorded logs which contain command detail and inferred results for either simulation/real vcf files. The default value is <em>distAngsdlog</em>. Such log files will always be provided.

* -method can be AmbiguityGT, NoAmbiguityGT, RandomSEQ, ConsensusSEQ, geno and nuc. The default value is <em>geno</em>.

* -model can be either JC or GTR. The default value is <em>JC</em>.

* -outglf specifies the directory of output genotype likelihoods files. If -outglf is not provided, genotype likelihoods files will not be saved. 

* -is2Dinfer determines whether genetic distance t and invariable site proportion p_inv are jointly inferred. The default value is 0, which means, by default, only genetic distance t will be inferred on the assumption that all sites are variable.
* -isex determines whether nucleotide transitions are excluded in the analysis. The default value is 0, which indicates nucleotide transitions are also included in the default analysis (It currently only work when -is2Dinfer is 0). 
* -p_inv gives the simulated true p_inv. It should be provided if a 2-d simulation and inference is conducted.

* -isthreading determines whether the EM algorithm is conducted parallelly. The default value is 0.

* -outbin determines whether the output files is saved as .bin or .txt.gz. The default value is 1.

* -outuchar determines whether store the output files in a more compressed unsinged char format. The default value is 0, which means the genotype likelihoods in vcf file will be, by default, stored in the double format.

* -inbin specifies the format of the input -inglf files (.bin or .txt.gz), the default value is 1. The default value is 1.

* -inuchar specifies whether the input -inglf files are stored in in a more compressed unsinged char format. The default value is 0.

* -numsites is the number of sites in simuation. The default number of sites is 1000000.

* -RD is the average read depth for simulation. The default value is 1.0.

* -e is the base calling error for simulation. The default value is 0.002.

* -tdiv is the true divergence time (genetic distance) for simulation. The default value is 1.0.

* -t1 -t2 are defined as manuscript. The default values of t1 and t2 are 0.4 and 0.25, correspondingly.

* -par gives comma-delimited 9 parameters (first 5 for symmetric matrix and last 4 for stationary distribution) for GTR model, and the last 4 as a distribution will be normalised automatically.

## Preparation and Command Examples
* Simulation
* distAngsd-geno: 
  * Input: vcf/bcf files of a pair of bam files.
  * Input Preparation:\
    bcftools mpileup -Ou -f ref.fa test1.sorted.bam test2.sorted.bam -o tmp.bcf;\
    bcftools filter -Ou -e INFO/INDEL!=0 tmp.bcf -o test.bcf;\
    bgzip test.bcf;\
    bcftools index test.bcf.gz;
  * Commands:
    * ./distAngsd -vcf test.bcf.gz [Options];
* distAngsd-geno: 
  * Input: genotype likelihood file [with either tab delimited txt.gz format or binary format]
  * Input Preparation:
    * txt.gz format: ./distAngsd -vcf test.bcf.gz -outglf test_nbin -outbin 0 [Options];
    * binary format: ./distAngsd -vcf test.bcf.gz -outglf test_bin -outbin 1 [Options];
  * Commands:
    * txt.gz format: ./distAngsd -inglf test_nbin.txt.gz -inbin 0 [Options];
    * binary format: ./distAngsd -inglf test_bin -inbin 1 [Options];
* distAngsd-nuc: 
  * Input: mpileup file of two bam files generated by samtools
  * Input Preparation:\
    samtools mpileup test1.sorted.bam test2.sorted.bam > test.txt;\
    bgzip test.txt;
  * Commands:\
    ./distAngsd -mpileup test.txt.gz -method nuc [Options];
    
## Cite
<b>Zhao, Lei#<b>; Nielsen, Rasmus; Korneliussen, Thorfinn S ; <a href="https://doi.org/10.1093/molbev/msac119"> distAngsd: Fast and accurate inference of genetic distances for Next Generation Sequencing data </a>, Molecular Biology and Evolution, 2022, 39(6), msac119, <a href="https://github.com/lz398/distAngsd"> Software Link </a>
