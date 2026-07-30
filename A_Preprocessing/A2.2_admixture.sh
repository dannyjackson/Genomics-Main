#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_file> -o <outprefix_for_files> -k <number_of_pops>

This script runs admixture from a VCF file (after plink conversion) for a series of inputted K values.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -v  Path to filtered vcf file
  -o  Output file prefix
  -k  Number of populations to test (defaults to 5). Will test all populations between 1 and your set value (ie: 1-5 if k=5)"
    exit 1
fi

NUM_K=5

# Parse command-line arguments
while getopts ":p:v:o:k:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCF=${OPTARG} ;;
        o) OUTPREFIX=${OPTARG} ;;
        k) NUM_K=${OPTARG} ;;
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
echo "Current script: A2.2_admixture.sh"

if [ -d "${OUTDIR}/analyses/admixture/" ];
        then
            echo "admixture directory already exists, moving on!"
        else
            echo "making admixture directory"
            mkdir -p "${OUTDIR}/analyses/admixture/"
fi

# First convert VCF to PLINK .bed format and convert chromids to ints
if [ -f "${OUTDIR}/analyses/admixture/${OUTPREFIX}.bed" ];
        then
            echo "PLINK .bed converted file for ${VCF} already exists, moving on!"
        else
            echo "Converting ${VCF} to PLINK .bed file"
            plink --vcf ${VCF} --allow-extra-chr --make-bed --out ${OUTDIR}/analyses/admixture/${OUTPREFIX}
            # Extra processing for admixture
            awk '{OFS="\t"; $1="0"; print}' ${OUTDIR}/analyses/admixture/${OUTPREFIX}.bim > ${OUTDIR}/analyses/admixture/${OUTPREFIX}_fixed.bim
            mv ${OUTDIR}/analyses/admixture/${OUTPREFIX}.bim ${OUTDIR}/analyses/admixture/${OUTPREFIX}_original.bim
            mv ${OUTDIR}/analyses/admixture/${OUTPREFIX}_fixed.bim ${OUTDIR}/analyses/admixture/${OUTPREFIX}.bim
fi

cd ${OUTDIR}/analyses/admixture/
echo ${NUM_K}
for k in $(seq ${NUM_K})
do 
    echo "Running admixture for K value: $k"
    admixture --cv ${OUTDIR}/analyses/admixture/${OUTPREFIX}.bed ${k} | tee admixlog_${OUTPREFIX}_${k}.out
    echo "Done running admixture for K $k"
done

echo "Done"


