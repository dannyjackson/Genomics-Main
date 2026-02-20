# msmc specific

source ~/.bashrc
micromamba activate msmc_env

source ../params_base.sh

# msmc input generation
k=150 # For snpability.sh
prefix=GCF_901933205 # Prefix of Reference Genome
MSMCTOOLS=${PROGDIR}/msmc-tools # directory with msmc-tools binaries
INPUTDIR=${OUTDIR}/datafiles/msmc/input # directory for msmc outputs
METHOD=samtools  # or another variant calling method

# msmc run params
NR_IND=5 # number of individuals in analysis
DATE=`date +%m%d%y`
RUN_NAME=msmc_${DATE}
P_PAR=1*2+25*1+1*2+1*3 
sex_chr=LG9 # name of sex chromosome to omit in analyses
nchr=`wc -l /xdisk/mcnew/finches/dannyjackson/finches/referencelists/SCAFFOLDS.txt | cut -f1 -d' '`
THREADS=28 # Number of CPUs to use in MSMC run. Should roughly match with the number of chromosomes you have
NUM_OPT=20 # Number of optimizations in MSMC run (20 is usually always enough)
INDEX="0,1,2,3,4,5,6,7,8,9,10,11" # This defaults to just the first 12 haplotype indices read per scaffold input.