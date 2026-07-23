#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -b <bamlist> -r <run_name> -c <chr_or_region_range>

This script calls variants for VCF file creation from BAM files for a single chromosome/region. 
Designed to be run in a slurm array for more efficient variant calling with larger sample sizes.
All outputs from array can be merged together using bcftools concat.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -b  Path to bam list file for analysis
  -r  Run name, required for providing a unique name to output files.
  -c  Either a Chromosome ID or region to analyze."
    exit 1
fi

# Parse command-line arguments
while getopts p:r:b:c: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        b) BAMLIST=${OPTARG};;
        c) CHR=${OPTARG};;
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
echo "Current script: A2.1_callvariants_array.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$REF" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi

bedtools makewindows -g ${REF}.fai -w 2000000 > ${OUTDIR}/referencelists/regions.bed

bcftools mpileup -Ou -f ${REF} -a FORMAT/AD,DP,INFO/AD,SP --bam-list ${BAMLIST} -r ${CHR} | bcftools call --threads ${THREADS} -mv -V indels -o ${OUTDIR}/datafiles/genotype_calls/"${RUNNAME}"_"${CHR}"_snps_multiallelic.vcf.gz
