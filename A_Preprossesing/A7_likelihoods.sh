#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -r <run_name> [-t]

This script computes average depth statistics of each sample from the output of A4_indelrealignment.sh.

Required arguments:
  -p  Path to the parameter file (e.g., params_preprocessing.sh).
  -r  Run name, required for providing a unique name to output files.

Optional argument:
  -t  Flag to exclude transition sites (-t activates -noTrans 1 in ANGSD).
"
    exit 1
}

# Check if no arguments were provided
if [ $# -lt 1 ]; then
    usage
fi

# Set default values for optional variables
TRANS_FLAG=""

# Parse command-line arguments
while getopts "p:r:t" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        r) RUNNAME=${OPTARG} ;;
        t) TRANS_FLAG="-noTrans 1" ;;
        *) usage ;;
    esac
done

# Ensure required arguments are provided
if [ -z "$PARAMS" ] || [ -z "$RUNNAME" ]; then
    echo "Error: Both -p (parameter file) and -r (run name) are required." >&2
    usage
fi

# Load parameters from the provided file
if [ ! -f "$PARAMS" ]; then
    echo "Error: Parameter file '$PARAMS' not found." >&2
    exit 1
fi

source "$PARAMS"

# Ensure required variables are set from the parameter file
REQUIRED_VARS=("OUTDIR" "THREADS" "REF" "FASTAS" "BAMUTILBAM" "SNPPVAL" "ANGSD" "PROJNAME")
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Missing required parameter '$var' in the parameter file." >&2
        exit 1
    fi


# Print script start information
echo -e "\n$(date)"
echo "Current script: A7_likelihoods.sh"

# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/geno_likelihoods"

# Generate genotype likelihoods
"${ANGSD}/angsd" -b "${OUTDIR}/referencelists/${PROJNAME}.bamlist.txt" \
  -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval "${SNPPVAL}" \
  -sites "${OUTDIR}/referencelists/${RUNNAME}.sites_headless.mafs" \
  -doBcf 1 -doGlf 2 -nThreads "${THREADS}" -out "${RUNNAME}.genolike" ${TRANS_FLAG}