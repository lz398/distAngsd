distAngsd is a software to simulate and infer phylogenetic relationships between two individuals.

make distAngsd

./distAngsd -o -method -model -glf -vcf -simrep -is2Dinfer -p_inv -isthreading -dobinary -numsites -isuchar -RD -errorrate -tdiv -t1 -t2 -par

-o specifies the recorded logs which contain command detail and inferred results for either simulation/real vcf files

-method can be either JC or GTR

-glf specifies the directory of output genotype likelihoods files

-vcf specifies the directory of input vcf files

-simrep specifies the number of simulation replicates

-is2Dinfer determines whether genetic distance t and invariable site proportion p_inv are jointly inferred

-p_inv gives the simulated true p_inv

-isthreading determines whether the EM algorithm is conducted parallelly

-dobinary determines whether the output files is saved as .bin or .txt.gz

-numsites is the number of sites in simuation

-isuchar determines whether store the real genotype likelihoods in a more compressed unsinged char format

-RD is the average read depth for simulation

-errorrate is the base calling error for simulation

-tdiv is the true divergence time (genetic distance) for simulation

-t1 -t2 are defined as manuscript

-par gives 9 parameters (first 5 for symmetric matrix and last 4 for stationary distribution) for GTR model, and the last 4 as a distribution will be normalised automatically.
