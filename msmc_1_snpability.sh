#!/bin/bash

# MSMC Pipeline Script
if [ $# -lt 1 ]
  then
    echo " This script generates a mask for a reference genome to determine which regions are "snpable".

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as params.sh"

  else
    while getopts p: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;

    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"

# Ensure required environment variables are set
if [ -z "$PROGDIR" ] || [ -z "$PROJHUB" ] || [ -z "$OUTDIR" ]; then
    echo "Error: PROGDIR, PROJHUB, or OUTDIR is not set." >&2
    exit 1
fi

snpable_script_path="${PROGDIR}/seqbility-20091110" # Directory with snpable scripts

# Check for required arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 -p <parameter_file> -g <reference_genome> -s <sample_codes_file>"
    exit 1
fi

while getopts ":p:g:s:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        g) GENOME=${OPTARG} ;;
        s) SAMPLES=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Ensure required arguments are provided
if [ -z "${PARAMS}" ] || [ -z "${GENOME}" ] || [ -z "${SAMPLES}" ]; then
    echo "Error: Missing required arguments." >&2
    exit 1
fi

# Load parameters from external file
source "${PARAMS}"

echo "Starting MSMC pipeline..."
date

# Clone necessary repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}"; exit 1; }
[ ! -d msmc2 ] && git clone https://github.com/stschiff/msmc2
[ ! -d msmc-tools ] && git clone https://github.com/stschiff/msmc-tools

# Install required dependencies
pip3 install --user whatshap

# Step 0: Create mappability mask
cd "${OUTDIR}" || { echo "Error: Could not change directory to ${OUTDIR}"; exit 1; }
export PATH=$PATH:$snpable_script_path

# Source additional parameters
source "${SCRIPTDIR}/msmc_params.sh"

mkdir -p "${OUTDIR}/snpable"
cd "${OUTDIR}/snpable" || { echo "Error: Could not change directory to ${OUTDIR}/snpable"; exit 1; }

echo "Starting extraction of overlapping ${k}-mer subsequences"

splitfa "${GENOME}" "${k}" | split -l 20000000
cat x* >> "${prefix}_split.${k}"

# Check if genome is indexed before proceeding
if [ ! -f "${GENOME}.bwt" ]; then
    echo "Indexing ${GENOME}"
    bwa index "${GENOME}"
else
    echo "${GENOME} already indexed."
fi

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to SAM file..."
bwa aln -t 8 -R 1000000 -O 3 -E 3 "${GENOME}" "${prefix}_split.${k}" > "${prefix}_split.${k}.sai"
bwa samse -f "${prefix}_split.${k}.sam" "${GENOME}" "${prefix}_split.${k}.sai" "${prefix}_split.${k}"

echo "Reads aligned, generating raw mask..."
gen_raw_mask.pl "${prefix}_split.${k}.sam" > "${prefix}_rawMask.${k}.fa"

echo "Generating final mask..."
gen_mask -l "${k}" -r 0.5 "${prefix}_rawMask.${k}.fa" > "${prefix}_mask.${k}.50.fa"

echo "Final mask saved as ${prefix}_mask.${k}.50.fa"


# Index reference genome
if ! command -v bwa &> /dev/null; then
    echo "Error: bwa is not installed." >&2
    exit 1
fi
bwa index "${GENOME}"

# Split genome into chunks
echo "Splitting genome..."
splitfa "${GENOME}" 150 | split -l 20000000
cat snpable/x* >> "${GENOME}_split.150"

# Align reads and generate mask
echo "Aligning reads to genome with BWA..."
bwa aln -t 8 -R 1000000 -O 3 -E 3 "${GENOME}" "${GENOME}_split.150" > "${GENOME}_split.150.sai"
bwa samse "${GENOME}" "${GENOME}_split.150.sai" "${GENOME}_split.150" > "${GENOME}_split.150.sam"

echo "Generating raw mask..."
gen_raw_mask.pl "${GENOME}_split.150.sam" > "${GENOME}_rawMask.150.fa"

echo "Generating final mask..."
gen_mask -l 150 -r 0.5 "${GENOME}_rawMask.150.fa" > "${GENOME}_mask.150.50.fa"

# Convert mask format
sed 's/>>/>/g' "${GENOME}_mask.150.50.fa" > "${GENOME}_revised_mask.150.50.fa"

# Create mappability mask
python2 "${PROGDIR}/msmc-tools/makeMappabilityMask.py"


