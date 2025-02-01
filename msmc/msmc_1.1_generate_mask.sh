#!/bin/bash

# MSMC Pipeline Script
# This script generates a mask for a reference genome to determine SNPable regions.

# Function to display usage
usage() {
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates a mask for a reference genome to determine SNPable regions."
    echo
    echo "Required Arguments:"
    echo "  -p  Path to the parameter file (e.g., params_msmc.sh from the GitHub repository)"
    exit 1
}

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    usage
fi

# Parse command-line arguments
while getopts "p:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        *) usage ;;
    esac
done

# Ensure parameter file is provided
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    usage
fi

# Load parameters from the provided file
if [ ! -f "${PARAMS}" ]; then
    echo "Error: Parameter file '${PARAMS}' not found!" >&2
    exit 1
fi

source "${PARAMS}"

# Ensure required environment variables are set
if [ -z "$PROGDIR" ] || [ -z "$PROJHUB" ] || [ -z "$OUTDIR" ]; then
    echo "Error: PROGDIR, PROJHUB, or OUTDIR is not set in the parameter file." >&2
    exit 1
fi

# Define path for Snpable scripts
SNPABLE_SCRIPT_PATH="${PROGDIR}/seqbility-20091110"

# Check if the Snpable directory exists
if [ ! -d "$SNPABLE_SCRIPT_PATH" ]; then
    echo "Error: Snpable script directory not found at ${SNPABLE_SCRIPT_PATH}." >&2
    exit 1
fi

# Add Snpable scripts to PATH
export PATH="$PATH:$SNPABLE_SCRIPT_PATH"

echo "Starting MSMC pipeline..."
date

# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d msmc2 ]; then
    git clone https://github.com/stschiff/msmc2 || { echo "Error: Failed to clone msmc2."; exit 1; }
fi

if [ ! -d msmc-tools ]; then
    git clone https://github.com/stschiff/msmc-tools || { echo "Error: Failed to clone msmc-tools."; exit 1; }
fi

echo "Environment setup completed."

# Step 0: Create mappability mask
cd "${OUTDIR}" || { echo "Error: Could not change directory to ${OUTDIR}."; exit 1; }

# Source additional MSMC parameters
if [ -f "${SCRIPTDIR}/msmc_params.sh" ]; then
    source "${SCRIPTDIR}/msmc_params.sh"
else
    echo "Error: MSMC parameter file '${SCRIPTDIR}/msmc_params.sh' not found." >&2
    exit 1
fi

mkdir -p "${OUTDIR}/snpable"
cd "${OUTDIR}/snpable" || { echo "Error: Could not change directory to ${OUTDIR}/snpable."; exit 1; }

echo "Extracting overlapping ${k}-mer subsequences..."
splitfa "${REF}" "${k}" | split -l 20000000
cat x* > "${prefix}_split.${k}"

# Check if REF is indexed before proceeding
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing ${REF} with BWA..."
    bwa index "${REF}" || { echo "Error: Failed to index ${REF}."; exit 1; }
else
    echo "Reference genome ${REF} is already indexed."
fi

# Align k-mer reads to the reference genome
echo "Aligning ${k}-mer reads to ${REF} with BWA..."
bwa aln -t 8 -R 1000000 -O 3 -E 3 "${REF}" "${prefix}_split.${k}" > "${prefix}_split.${k}.sai" || { echo "Error: BWA alignment failed."; exit 1; }

bwa samse -f "${prefix}_split.${k}.sam" "${REF}" "${prefix}_split.${k}.sai" "${prefix}_split.${k}" || { echo "Error: BWA SAM conversion failed."; exit 1; }

echo "Reads aligned. Generating raw mask..."
gen_raw_mask.pl "${prefix}_split.${k}.sam" > "${prefix}_rawMask.${k}.fa" || { echo "Error: Failed to generate raw mask."; exit 1; }

echo "Generating final mask..."
gen_mask -l "${k}" -r 0.5 "${prefix}_rawMask.${k}.fa" > "${prefix}_mask.${k}.50.fa" || { echo "Error: Failed to generate final mask."; exit 1; }

echo "Final mask saved as ${prefix}_mask.${k}.50.fa"
echo "MSMC pipeline completed successfully!"
date
