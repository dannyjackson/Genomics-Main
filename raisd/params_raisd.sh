# params_raisd.sh

source ${SCRIPTDIR}/Genomics-Main/params_base.sh

CUTOFF=threshold_for_top_genes_or_snps_eg_0.01

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1=#XXXXXX
COLOR2=#XXXXXX

REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna # path to reference genome
SCAFFOLD_LIST=/xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt

# source the setup file for raisd                                           
source ${SCRIPTDIR}/Genomics-Main/raisd/setup_raisd.sh