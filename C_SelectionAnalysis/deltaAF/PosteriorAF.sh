#!/bin/sh

# FST computation script using ANGSD
# Computes Fst between two genome groups using genotype likelihoods

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $(basename $0) -p <param_file> [-w <window_size>] [-s <step_size>]"
    echo "\nThis script computes Fst between two genome groups using ANGSD."
    echo "It requires SAF files as input and produces genome-wide Fst estimates,"
    echo "as well as sliding window Fst statistics.\n"
    echo "Required Argument:"
    echo "  -p   Path to parameter file (see example params_fst.sh in the repository)\n"
    echo "  -s   Name of species"
    echo "  -q   Name of population"
    exit 1
fi

# Parse arguments
while getopts "p:s:q:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        s) SP={OPTARG} ;;
        q) POP={OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

# Ensure required parameter file is provided
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

# Set defaults if not provided
WIN=${WIN:-10000}
STEP=${STEP:-10000}

# Read chromosome file 
CHROM=$(cat "$CHR_FILE")


# Print script info
echo "\nRunning FST computation script"


# Re-run with posterior MAFs using the SFS as prior
mkdir -p ${OUTDIR}/analyses/deltaAF/${SP}/

${ANGSD}/angsd -b /xdisk/mcnew/finches/dannyjackson/finches/referencelists/${SP}${POP}bams.txt \
  -ref /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna \
  -out ${OUTDIR}/analyses/deltaAF/${SP}/${SP}_${POP}_postmafs \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
  -minMapQ 30 -minQ 20 \
  -GL 1 \
  -doMaf 1 -doMajorMinor 1 -doPost 1 \
  -pest ${OUTDIR}/datafiles/safs/${SP}${POP}.sfs \
  -minInd 3