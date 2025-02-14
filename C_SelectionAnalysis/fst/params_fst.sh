# params_fst
source /path/to/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=threshold_for_top_genes_or_snps_eg_0.01

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1=#XXXXXX
COLOR2=#XXXXXX

# define the names of the two populations that will be compared
POP1=population1name
POP2=population2name

# source the setup file for fst
source ${SCRIPTDIR}/Genomics-Main/fst/setup_fst.sh
