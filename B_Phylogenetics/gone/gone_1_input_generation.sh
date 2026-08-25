n#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

Use VCFtools and PLINK 1.9 to generate some basic input files for GONE2
This step often requires trial-and-error with subsetting chromosomes/samples/snps and various filtering of missing data and MAF. Customize as needed
You may consider looking into the following filtering steps if gone throws issues about missing data or memory requirements
    VCFTOOLS: --max-missing
    PLINK: --thin-count / --maf

Required argument:
  -p  Path to the main parameter file (e.g., params_base.sh in the GitHub repository).
  -v  Path to raw VCF file
  -s  Path to file containing scaffold subset to include
  -o  Output name to assign to files"
    exit 1
fi

# Parse command-line arguments
while getopts p:v:s:o: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCF=${OPTARG};;
        s) SCAFFOLD_LIST=${OPTARG};;
        o) OUTNAME=${OPTARG};;
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
echo "Current script: gone_1_input_generation.sh"

CHR_SUBSET=$(for name in $(cat $SCAFFOLD_LIST); do echo --chr $name; done)

RESULT_DIR=${OUTDIR}/datafiles/gone2_inputs/${OUTNAME}
if [ ! -d "$RESULT_DIR" ]; then
  echo "Directory for gone2 input for ${OUTNAME} does not exist. Creating it now..."
  mkdir -p "$RESULT_DIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone2 input for ${OUTNAME} already exists. Moving on..."
fi

# Create a filtered VCF to only include a specific subset of chromosomes
vcftools $CHR_SUBSET --vcf $VCF --recode --recode-INFO-all --out ${RESULT_DIR}/${OUTNAME}
mv ${RESULT_DIR}/${OUTNAME}.recode.vcf ${RESULT_DIR}/${OUTNAME}.vcf 

# Convert VCF to plink formats
plink --vcf ${RESULT_DIR}/${OUTNAME}.vcf --allow-extra-chr --recode --out ${RESULT_DIR}/${POPNAME}
# NOTE: PLINK .map file will not generate with required cM distances for GONE2 analysis.