#!/bin/bash
#run locally like
#make test BAMDIR=smallBam/


PRG=../distAngsd
echo "--------------------"
echo "Using PRG: '${PRG}'"
echo "--------------------"
#JC model
echo "Switch to JC model"
#VCF
echo "Testing VCF/BCF Input disAngsd-geno"
echo "Binary output"
$PRG -vcf test.bcf.gz;
md5sum distAngsdlog > JC/vcfb.md5;
$PRG -vcf test.bcf.gz -isex 1;
md5sum distAngsdlog > JC/vcfbex.md5;
echo "Text output"
$PRG -vcf test.bcf.gz -outbin 0;
md5sum distAngsdlog > JC/vcfnb.md5;
$PRG -vcf test.bcf.gz -outbin 0 -isex 1;
md5sum distAngsdlog > JC/vcfnbex.md5;
echo "Binary output + threading"
$PRG -vcf test.bcf.gz -isthreading 1;
md5sum distAngsdlog > JC/vcfbth.md5;
$PRG -vcf test.bcf.gz -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/vcfbthex.md5;
#md5sum -c JC/vcfbth.md5 || exit 1;
echo "Text output + threading"
$PRG -vcf test.bcf.gz -outbin 0 -isthreading 1;
md5sum distAngsdlog > JC/vcfnbth.md5;
$PRG -vcf test.bcf.gz -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/vcfnbthex.md5;
#md5sum -c JC/vcfnbth.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -inglf test_bin;
md5sum distAngsdlog > JC/glfbb.md5;
$PRG -inglf test_bin -isex 1;
md5sum distAngsdlog > JC/glfbbex.md5;
#md5sum -c JC/glfbb.md5 || exit 1;
echo "Binary input text output"
$PRG -inglf test_bin -outbin 0;
md5sum distAngsdlog > JC/glfbnb.md5;
$PRG -inglf test_bin -outbin 0 -isex 1;
md5sum distAngsdlog > JC/glfbnbex.md5;
#md5sum -c JC/glfbnb.md5 || exit 1;
echo "Text input binary output"
$PRG -inglf test_nbin.txt.gz -inbin 0;
md5sum distAngsdlog > JC/glfnbb.md5;
$PRG -inglf test_nbin.txt.gz -inbin 0 -isex 1;
md5sum distAngsdlog > JC/glfnbbex.md5;
#md5sum -c JC/glfnbb.md5 || exit 1;
echo "Text input text output"
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0;
md5sum distAngsdlog > JC/glfnbnb.md5;
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isex 1;
md5sum distAngsdlog > JC/glfnbnbex.md5;
#md5sum -c JC/glfnbnb.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -inglf test_bin -isthreading 1;
md5sum distAngsdlog > JC/glfbbth.md5;
$PRG -inglf test_bin -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/glfbbthex.md5;
#md5sum -c JC/glfbbth.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -inglf test_bin -outbin 0 -isthreading 1;
md5sum distAngsdlog > JC/glfbnbth.md5;
$PRG -inglf test_bin -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/glfbnbthex.md5;
#md5sum -c JC/glfbnbth.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -inglf test_nbin.txt.gz -inbin 0 -isthreading 1;
md5sum distAngsdlog > JC/glfnbbth.md5;
$PRG -inglf test_nbin.txt.gz -inbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/glfnbbthex.md5;
#md5sum -c JC/glfnbbth.md5 || exit 1;
echo "Text input text output + threading"
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isthreading 1;
md5sum distAngsdlog > JC/glfnbnbth.md5;
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/glfnbnbthex.md5;
#md5sum -c JC/glfnbnbth.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -mpileup test.txt.gz;
md5sum distAngsdlog > JC/mpileupb.md5;
$PRG -method nuc -mpileup test.txt.gz -isex 1;
md5sum distAngsdlog > JC/mpileupbex.md5;
#md5sum -c JC/mpileupb.md5 || exit 1;
echo "Text output"
$PRG -method nuc -mpileup test.txt.gz -outbin 0;
md5sum distAngsdlog > JC/mpileupnb.md5;
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isex 1;
md5sum distAngsdlog > JC/mpileupnbex.md5;
#md5sum -c JC/mpileupnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -mpileup test.txt.gz -isthreading 1;
md5sum distAngsdlog > JC/mpileupbth.md5;
$PRG -method nuc -mpileup test.txt.gz -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/mpileupbthex.md5;
#md5sum -c JC/mpileupbth.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isthreading 1;
md5sum distAngsdlog > JC/mpileupnbth.md5;
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > JC/mpileupnbthex.md5;
#md5sum -c JC/mpileupnbth.md5 || exit 1;

#GTR model
echo "Switch to GTR model"
#VCF
echo "Testing VCF/BCF Input disAngsd-geno"
echo "Binary output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum distAngsdlog > GTR/vcfb.md5;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum distAngsdlog > GTR/vcfbex.md5;
#md5sum -c GTR/vcfb.md5 || exit 1;
echo "Text output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum distAngsdlog > GTR/vcfnb.md5;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum distAngsdlog > GTR/vcfnbex.md5;
#md5sum -c GTR/vcfnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum distAngsdlog > GTR/vcfbth.md5;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/vcfbthex.md5;
#md5sum -c GTR/vcfbth.md5 || exit 1;
echo "Text output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR/vcfnbth.md5;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/vcfnbthex.md5;
#md5sum -c GTR/vcfnbth.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum distAngsdlog > GTR/glfbb.md5;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum distAngsdlog > GTR/glfbbex.md5;
#md5sum -c GTR/glfbb.md5 || exit 1;
echo "Binary input text output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum distAngsdlog > GTR/glfbnb.md5;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum distAngsdlog > GTR/glfbnbex.md5;
#md5sum -c GTR/glfbnb.md5 || exit 1;
echo "Text input binary output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum distAngsdlog > GTR/glfnbb.md5;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum distAngsdlog > GTR/glfnbbex.md5;
#md5sum -c GTR/glfnbb.md5 || exit 1;
echo "Text input text output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum distAngsdlog > GTR/glfnbnb.md5;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum distAngsdlog > GTR/glfnbnbex.md5;
#md5sum -c GTR/glfnbnb.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum distAngsdlog > GTR/glfbbth.md5;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/glfbbthex.md5;
#md5sum -c GTR/glfbbth.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR/glfbnbth.md5;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/glfbnbthex.md5;
#md5sum -c GTR/glfbnbth.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum distAngsdlog > GTR/glfnbbth.md5;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/glfnbbthex.md5;
#md5sum -c GTR/glfnbbth.md5 || exit 1;
echo "Text input text output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR/glfnbnbth.md5;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/glfnbnbthex.md5;
#md5sum -c GTR/glfnbnbth.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum distAngsdlog > GTR/mpileupb.md5;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum distAngsdlog > GTR/mpileupbex.md5;
#md5sum -c GTR/mpileupb.md5 || exit 1;
echo "Text output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum distAngsdlog > GTR/mpileupnb.md5;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum distAngsdlog > GTR/mpileupnbex.md5;
#md5sum -c GTR/mpileupnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum distAngsdlog > GTR/mpileupbth.md5;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/mpileupbthex.md5;
#md5sum -c GTR/mpileupbth.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR/mpileupnbth.md5;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum distAngsdlog > GTR/mpileupnbthex.md5;
#md5sum -c GTR/mpileupnbth.md5 || exit 1;

#GTR+I model
echo "Switch to GTR+I model"
#VCF
echo "Testing VCF/BCF Input disAngsd-geno"
echo "Binary output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum distAngsdlog > GTR2D/vcfb.md5;
#md5sum -c GTR2D/vcfb.md5 || exit 1;
echo "Text output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum distAngsdlog > GTR2D/vcfnb.md5;
#md5sum -c GTR2D/vcfnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum distAngsdlog > GTR2D/vcfbth.md5;
#md5sum -c GTR2D/vcfbth.md5 || exit 1;
echo "Text output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR2D/vcfnbth.md5;
#md5sum -c GTR2D/vcfnbth.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum distAngsdlog > GTR2D/glfbb.md5;
#md5sum -c GTR2D/glfbb.md5 || exit 1;
echo "Binary input text output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum distAngsdlog > GTR2D/glfbnb.md5;
#md5sum -c GTR2D/glfbnb.md5 || exit 1;
echo "Text input binary output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum distAngsdlog > GTR2D/glfnbb.md5;
#md5sum -c GTR2D/glfnbb.md5 || exit 1;
echo "Text input text output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum distAngsdlog > GTR2D/glfnbnb.md5;
#md5sum -c GTR2D/glfnbnb.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum distAngsdlog > GTR2D/glfbbth.md5;
#md5sum -c GTR2D/glfbbth.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR2D/glfbnbth.md5;
#md5sum -c GTR2D/glfbnbth.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum distAngsdlog > GTR2D/glfnbbth.md5;
#md5sum -c GTR2D/glfnbbth.md5 || exit 1;
echo "Text input text output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR2D/glfnbnbth.md5;
#md5sum -c GTR2D/glfnbnbth.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum distAngsdlog > GTR2D/mpileupb.md5;
#md5sum -c GTR2D/mpileupb.md5 || exit 1;
echo "Text output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum distAngsdlog > GTR2D/mpileupnb.md5;
#md5sum -c GTR2D/mpileupnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum distAngsdlog > GTR2D/mpileupbth.md5;
#md5sum -c GTR2D/mpileupbth.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum distAngsdlog > GTR2D/mpileupnbth.md5;
#md5sum -c GTR2D/mpileupnbth.md5 || exit 1;
