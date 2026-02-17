#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This uses VCFtools to generate input files for LDHelmet

Required argument:
  -p  Path to the parameter file (e.g., params_ldhelmet.sh in the GitHub repository).
  -v  VCF file path. Must phased. Assuming bgzipped
  -c  Path to file containing chromosome names to proceess, one per line."
    exit 1
fi

# Parse command-line arguments
while getopts p:v:c: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCF=${OPTARG};;
        c) CHR_LST=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

printf "\n\n\n\n"
date
echo "Current script: ldhelmet_input_generation.sh"

module load vcftools

OUTPREFIX=$(basename "$VCF" .vcf.gz)

if [ ! -d "${INPUT_DIR}/${OUTPREFIX}" ]; then
        echo "${OUTPREFIX} input directory already exists."
    else
        echo "${OUTPREFIX} input directory does not exist. Creating it now..."
        mkdir -p "${INPUT_DIR}/${OUTPREFIX}"
    fi

for CHR in $(cat ${CHR_LST}); do
    echo "Processing chromosome: ${CHR}"

    vcftools --gzvcf ${VCF} --chr ${CHR} --ldhelmet --out ${INPUT_DIR}/${OUTPREFIX}/${CHR}
done