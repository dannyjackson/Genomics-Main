# params_raisd.sh

source /path/to/params_base.sh

CUTOFF=threshold_for_top_genes_or_snps_eg_0.01

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1=#XXXXXX
COLOR2=#XXXXXX

# define the name of the  population that will be analyzed
POP1=population1name


# source the setup file for raisd                                           
source ${SCRIPTDIR}/Genomics-Main/raisd/setup_raisd.sh