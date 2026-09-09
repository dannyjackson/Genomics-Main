source ../params_base.sh
source ${PROGDIR}/.flexsweep_venv/bin/activate

# General
POPNAME=pop1
THREADS=24

# Simulations
NUM_HAPS=XX # Number of haplotypes (2x sample size)
DEMES=/path/to/demes/yaml/file
SIMULATIONS=250000

# VCF feature vectors
VCFDIR=/path/to/vcf/directory
REC_MAP=/path/to/recombination/map/file