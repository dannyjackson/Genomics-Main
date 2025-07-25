# msmc specific
# msmc input generation
k=150 # For snpability.sh
prefix=GCF_901933205
MSMCTOOLS=${PROGDIR}/msmc-tools # directory with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS # add directory with msmc-tools binaries to path
METHOD=samtools  # or another variant calling method

# msmc run params
NR_IND=1 # number of individuals in analysis
POP_OR_IND=SRR2917338 # name of individual or population being analyzed for script 3
DATE=`date +%m%d%y`
RUN_NAME=msmc_${DATE}
P_PAR=1*2+25*1+1*2+1*3 
sex_chr=LG9 # name of sex chromosome to omit in analyses
nchr=`wc -l /xdisk/mcnew/finches/dannyjackson/finches/referencelists/SCAFFOLDS.txt | cut -f1 -d' '`