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
md5sum -c JC/vcfb.md5 || exit 1;
$PRG -vcf test.bcf.gz -isex 1;
md5sum -c JC/vcfbex.md5 || exit 1;
echo "Text output"
$PRG -vcf test.bcf.gz -outbin 0;
md5sum -c JC/vcfnb.md5 || exit 1;
$PRG -vcf test.bcf.gz -outbin 0 -isex 1;
md5sum -c JC/vcfnbex.md5 || exit 1;
echo "Binary output + threading"
$PRG -vcf test.bcf.gz -isthreading 1;
md5sum -c JC/vcfbth.md5 || exit 1;
$PRG -vcf test.bcf.gz -isthreading 1 -isex 1;
md5sum -c JC/vcfbthex.md5 || exit 1;
echo "Text output + threading"
$PRG -vcf test.bcf.gz -outbin 0 -isthreading 1;
md5sum -c JC/vcfnbth.md5 || exit 1;
$PRG -vcf test.bcf.gz -outbin 0 -isthreading 1 -isex 1;
md5sum -c JC/vcfnbthex.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -inglf test_bin;
md5sum -c JC/glfbb.md5 || exit 1;
$PRG -inglf test_bin -isex 1;
md5sum -c JC/glfbbex.md5 || exit 1;
echo "Binary input text output"
$PRG -inglf test_bin -outbin 0;
md5sum -c JC/glfbnb.md5 || exit 1;
$PRG -inglf test_bin -outbin 0 -isex 1;
md5sum -c JC/glfbnbex.md5 || exit 1;
echo "Text input binary output"
$PRG -inglf test_nbin.txt.gz -inbin 0;
md5sum -c JC/glfnbb.md5 || exit 1;
$PRG -inglf test_nbin.txt.gz -inbin 0 -isex 1;
md5sum -c JC/glfnbbex.md5 || exit 1;
echo "Text input text output"
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0;
md5sum -c JC/glfnbnb.md5 || exit 1;
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isex 1;
md5sum -c JC/glfnbnbex.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -inglf test_bin -isthreading 1;
md5sum -c JC/glfbbth.md5 || exit 1;
$PRG -inglf test_bin -isthreading 1 -isex 1;
md5sum -c JC/glfbbthex.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -inglf test_bin -outbin 0 -isthreading 1;
md5sum -c JC/glfbnbth.md5 || exit 1;
$PRG -inglf test_bin -outbin 0 -isthreading 1 -isex 1;
md5sum -c JC/glfbnbthex.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -inglf test_nbin.txt.gz -inbin 0 -isthreading 1;
md5sum -c JC/glfnbbth.md5 || exit 1;
$PRG -inglf test_nbin.txt.gz -inbin 0 -isthreading 1 -isex 1;
md5sum -c JC/glfnbbthex.md5 || exit 1;
echo "Text input text output + threading"
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isthreading 1;
md5sum -c JC/glfnbnbth.md5 || exit 1;
$PRG -inglf test_nbin.txt.gz -inbin 0 -outbin 0 -isthreading 1 -isex 1;
md5sum -c JC/glfnbnbthex.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -mpileup test.txt.gz;
md5sum -c JC/mpileupb.md5 || exit 1;
$PRG -method nuc -mpileup test.txt.gz -isex 1;
md5sum -c JC/mpileupbex.md5 || exit 1;
echo "Text output"
$PRG -method nuc -mpileup test.txt.gz -outbin 0;
md5sum -c JC/mpileupnb.md5 || exit 1;
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isex 1;
md5sum -c JC/mpileupnbex.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -mpileup test.txt.gz -isthreading 1;
md5sum -c JC/mpileupbth.md5 || exit 1;
$PRG -method nuc -mpileup test.txt.gz -isthreading 1 -isex 1;
md5sum -c JC/mpileupbthex.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isthreading 1;
md5sum -c JC/mpileupnbth.md5 || exit 1;
$PRG -method nuc -mpileup test.txt.gz -outbin 0 -isthreading 1 -isex 1;
md5sum -c JC/mpileupnbthex.md5 || exit 1;

#GTR model
echo "Switch to GTR model"
#VCF
echo "Testing VCF/BCF Input disAngsd-geno"
echo "Binary output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum -c GTR/vcfb.md5 || exit 1;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum -c GTR/vcfbex.md5 || exit 1;
echo "Text output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum -c GTR/vcfnb.md5 || exit 1;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum -c GTR/vcfnbex.md5 || exit 1;
echo "Binary output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum -c GTR/vcfbth.md5 || exit 1;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum -c GTR/vcfbthex.md5 || exit 1;
echo "Text output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum -c GTR/vcfnbth.md5 || exit 1;
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum -c GTR/vcfnbthex.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum -c GTR/glfbb.md5 || exit 1;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum -c GTR/glfbbex.md5 || exit 1;
echo "Binary input text output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum -c GTR/glfbnb.md5 || exit 1;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum -c GTR/glfbnbex.md5 || exit 1;
echo "Text input binary output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum -c GTR/glfnbb.md5 || exit 1;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum -c GTR/glfnbbex.md5 || exit 1;
echo "Text input text output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum -c GTR/glfnbnb.md5 || exit 1;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum -c GTR/glfnbnbex.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum -c GTR/glfbbth.md5 || exit 1;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum -c GTR/glfbbthex.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum -c GTR/glfbnbth.md5 || exit 1;
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum -c GTR/glfbnbthex.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum -c GTR/glfnbbth.md5 || exit 1;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum -c GTR/glfnbbthex.md5 || exit 1;
echo "Text input text output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum -c GTR/glfnbnbth.md5 || exit 1;
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum -c GTR/glfnbnbthex.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946;
md5sum -c GTR/mpileupb.md5 || exit 1;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isex 1;
md5sum -c GTR/mpileupbex.md5 || exit 1;
echo "Text output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0;
md5sum -c GTR/mpileupnb.md5 || exit 1;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isex 1;
md5sum -c GTR/mpileupnbex.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1;
md5sum -c GTR/mpileupbth.md5 || exit 1;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -isthreading 1 -isex 1;
md5sum -c GTR/mpileupbthex.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1;
md5sum -c GTR/mpileupnbth.md5 || exit 1;
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -outbin 0 -isthreading 1 -isex 1;
md5sum -c GTR/mpileupnbthex.md5 || exit 1;

#GTR+I model
echo "Switch to GTR+I model"
#VCF
echo "Testing VCF/BCF Input disAngsd-geno"
echo "Binary output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum -c GTR2D/vcfb.md5 || exit 1;
echo "Text output"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum -c GTR2D/vcfnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum -c GTR2D/vcfbth.md5 || exit 1;
echo "Text output + threading"
$PRG -model GTR -vcf test.bcf.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum -c GTR2D/vcfnbth.md5 || exit 1;

#GLF
echo "Testing Genotype likelihood Input disAngsd-geno"
echo "Binary input binary output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum -c GTR2D/glfbb.md5 || exit 1;
echo "Binary input text output"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum -c GTR2D/glfbnb.md5 || exit 1;
echo "Text input binary output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum -c GTR2D/glfnbb.md5 || exit 1;
echo "Text input text output"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum -c GTR2D/glfnbnb.md5 || exit 1;
echo "Binary input binary output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum -c GTR2D/glfbbth.md5 || exit 1;
echo "Binary input text output + threading"
$PRG -model GTR -inglf test_bin -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum -c GTR2D/glfbnbth.md5 || exit 1;
echo "Text input binary output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum -c GTR2D/glfnbbth.md5 || exit 1;
echo "Text input text output + threading"
$PRG -model GTR -inglf test_nbin.txt.gz -inbin 0 -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum -c GTR2D/glfnbnbth.md5 || exit 1;

#MPILEUP
echo "Testing mpileup Input disAngsd-nuc"
echo "Binary output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1;
md5sum -c GTR2D/mpileupb.md5 || exit 1;
echo "Text output"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0;
md5sum -c GTR2D/mpileupnb.md5 || exit 1;
echo "Binary output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -isthreading 1;
md5sum -c GTR2D/mpileupbth.md5 || exit 1;
echo "Text output + threading"
$PRG -method nuc -model GTR -mpileup test.txt.gz -par 2.0431,0.0821,0,0.067,0,0.2184,0.2606,0.3265,0.1946 -is2Dinfer 1 -outbin 0 -isthreading 1;
md5sum -c GTR2D/mpileupnbth.md5 || exit 1;
#PRG=""
#BAMDIR=""
#
#if [ $# -eq 0 ]
#then
#    exit 1;
#fi
#
#if [ $# -eq 1 ]
#then
#    PRG=$1
#fi
#
#if [ $# -eq 2 ]
#then
#    PRG=$1
#    BAMDIR=$2
#fi
#
#
#echo "--------------------"
#echo "Using PRG: '${PRG}' and BAMDIR: '${BAMDIR}'"
#echo "--------------------"
#
#WDIR=`dirname $PRG`
#
#RVAL=0
#
#
#if [[ ! -z "$BAMDIR" ]]; then
#echo "Testing -sites"
#time ./testVcf.sh $WDIR/angsd $BAMDIR
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with -sites exit code: $?"
##    cat ./testFilterSites.sh.log
#    RVAL=1
#fi
#fi
#
#if [[ ! -z "$BAMDIR" ]]; then
#echo "Testing vcfreading"
#time ./testVcf.sh $WDIR/angsd ${BAMDIR}/small2.bcf
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with -vcf-pl exit code: $?"
#    cat ./testVcf.sh.log
#    RVAL=1
#fi
#fi
#
#echo "Testing neutrality test statistics"
#time ./testTaj.sh $WDIR
#if [ ! $? -eq 0 ] ;then
#    echo "Problem with neutrality test statistics exit code: $?"
#    cat ./testTaj.sh.log
#    RVAL=1
#fi
#
#echo "Testing neutrality test statistics (haploid)"
#time ./testHapSFS.sh $WDIR
#if [ ! $? -eq 0 ] ;then
#    echo "Problem with neutrality test statistics exit code: $?"
#    cat ./testHapSFS.sh.log
#    RVAL=1
#fi
#
#echo "Testing fst using msms"
#time ./testFst.sh $WDIR
#if [ ! $? -eq 0 ] ;then
#    echo "Problem with Fst test statistics exit code: $?"
#    cat ./testFst.sh.log
#    RVAL=1
#fi
#
#echo "Testing fst_folded using msms"
#time ./testFst_folded.sh $WDIR
#if [ ! $? -eq 0 ] ;then
#    echo "Problem with Fst_folded test statistics exit code: $?"
#    cat ./testFst_folded.sh.log
#    RVAL=1
#fi
#
#
#
#echo "Testing SFS"
#time ./testSFS.sh $WDIR
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with SFS exit code: $?"
#    cat ./testSFS.sh.log
#    RVAL=1
#fi
#
#if [[ ! -z "$BAMDIR" ]]; then
#echo "Testing basic mpileup"
#time ./testBam.sh $WDIR/angsd $BAMDIR
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with basic pileup exit code: $?"
#    cat ./testBam.sh.log
#    RVAL=1
#fi
#fi
#
#echo "Testing association"
#time ./testDoAsso2456.sh $WDIR
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with association exit code: $?"
#    cat ./testDoAsso2456.sh.log
#    RVAL=1
#fi
#
#if [[ ! -z "$BAMDIR" ]]; then
#echo "Testing haplocall"
#time ./testHaploCall.sh $WDIR/angsd
#if [ ! $? -eq 0  ]   ;then
#    echo "Problem with haplocall exit code: $?"
#    cat ./testHaploCall.sh.log
#    RVAL=1
#fi
#fi
#
#
#exit ${RVAL}
#
#if false; then
#
#    echo testAsso6.sh $PRG
#    ./testAsso6.sh $PRG
#
#    #echo ./testErr.sh $PRG
#    #./testErr.sh $PRG
#
#    echo ./testGL6.sh $PRG
#    ./testGL6.sh $PRG
#
#    echo ./testMisc9.sh $PRG
#    ./testMisc9.sh $PRG
#
#    echo ./testAbba.sh $PRG
#    ./testAbba.sh $PRG
#
#    echo ./testFasta.sh
#    ./testFasta.sh $PRG
#
#    echo ./testBaq.sh $PRG
#    ./testBaq3.sh $PRG
#
#    echo "Netaccess is now deprecated"
#    #echo ./testNetAccess.sh $PRG
#    #./testNetAccess.sh $PRG
#fi
