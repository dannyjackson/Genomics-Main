# msmc specific

source ~/.bashrc
micromamba activate msmc_env

source ../../params_base.sh

# input generation
POPNAME=pop1
INDFILE=${OUTDIR}/referencelists/pop1_sampleids.txt

# msmc analysis
NR_IND=5 # number of individuals in analysis
P_PAR=1*2+25*1+1*2+1*3 
THREADS=28 # Number of CPUs. Should roughly match with the number of chromosomes you have
NUM_OPT=20 # Number of optimizations (20 is usually enough)
INDEX="0,1,2,3,4,5,6,7,8,9,10,11" # This defaults to just the first 12 haplotype indices read per scaffold input.