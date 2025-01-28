module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2

MSMCTOOLS=~/programs/msmc-tools # directory with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS

OUTDIR=/xdisk/mcnew/dannyjackson/cardinals_dfinch/    # main directory for output files
ANGSD=~/programs/angsd/     # path to directory with angsd executables

scriptdir=~/programs/Genomics-Main/
PATH=$PATH:$scriptdir # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path

# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/fst
mkdir -p ${OUTDIR}/analyses/genelist
mkdir -p ${OUTDIR}/datafiles/safs
mkdir -p ${OUTDIR}/datafiles/mls/
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
mkdir -p ${OUTDIR}/analyses/genelist/${WIN}


# for fst_1.sh

# define the names of the two populations that will be compared
POP1=nocaurban
POP2=nocarural

# characters at the start of a chromosome number (excluding scaffolds)
# note that this script also assumes chromosomes will end in .1 i.e. NC_012345.1, and will remove the .1 from files where it will disrupt plotting etc

CHRLEAD=NC_0
# path to reference genome
REF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
# path to gff file
GFF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCA_013397215.1/genomic.gff 

