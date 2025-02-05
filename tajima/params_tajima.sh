# params_tajima
source /path/to/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=threshold_for_top_genes_or_snps_eg_0.01

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1=#XXXXXX
COLOR2=#XXXXXX

# define the name of the population that will be analyzed
POP=populationname

source ${SCRIPTDIR}/Genomics-Main/tajima/setup_tajima.sh