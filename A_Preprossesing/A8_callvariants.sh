#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script computes average depth statistics of each bam file in a directory.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -b  Path to bam directory for analysis (if you are following the full pipeline, this will be OUTDIR/datafiles/indelrealignment/)
  -r  Run name, required for providing a unique name to output files."
    exit 1
fi

# Parse command-line arguments
while getopts p:b:r: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        b) BAMDIR=${OPTARG};;
        r) RUNNAME=${OPTARG};;
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
echo "Current script: A2_ClipOverlap.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$REF" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/genotype_calls/"


bcftools mpileup -Ou -f ${REF} -a FORMAT/AD,DP,INFO/AD,SP "${BAMDIR}"*.bam | bcftools call -mv -V indels > ${OUTDIR}/datafiles/genotype_calls/"$RUNNAME"_snps_multiallelic.vcf

