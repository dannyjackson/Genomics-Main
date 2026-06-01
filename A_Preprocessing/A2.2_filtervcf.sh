#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script filters the VCF generated in A2.1_callvariants.sh"
    echo
    echo "Required Arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh from the GitHub repository)"
    echo "  -i  Prefix of vcf file"
    echo "  -f  optional boolean flags to fill extra info tags (allele frequency, count, and number) Required for SAPPHIRE Rephasing (A2.7)"
    exit 1
}

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    usage
fi

FILLTAGS=false

# Parse command-line arguments
while getopts "p:i:f" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) ID=${OPTARG} ;;
        f) FILLTAGS=true ;;
        *) usage ;;
    esac
done

source ${PARAMS}


# Filter VCF based on quality
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf" ]
        then
            echo "qualitysort vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by quality greater than 100"
            bcftools view -i 'QUAL>100' "${OUTDIR}/datafiles/genotype_calls/${ID}_snps_multiallelic.vcf.gz" > "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf"

fi


# Filter VCF based on depth, remove indels
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_depthfiltered.vcf" ]
        then
            echo "depthfilter vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by depth, remove indels"
            vcftools --vcf "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf" --min-meanDP ${MINDEPTH} --max-meanDP ${MAXDEPTH} --remove-indels --recode --out "${OUTDIR}/datafiles/genotype_calls/${ID}_depthfiltered"
fi


# Further filtering using PLINK
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_plinkfiltered.vcf" ]
        then
            echo "plink filtered vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by geno 0.2, minimum ind 0.2, maf 0.01, remove indels"
            plink --vcf "${OUTDIR}/datafiles/genotype_calls/${ID}_depthfiltered.recode.vcf" \
            --allow-extra-chr --snps-only 'just-acgt' \
            --geno 0.02 --mind 0.2 --maf 0.01 \
            --recode vcf-iid --out "${OUTDIR}/datafiles/genotype_calls/${ID}_plinkfiltered"
fi

# Filling INFO tags for AD, AC, and AN and generating BCF
if [ "$FILLTAGS" = true ]; then
    echo "Filling INFO tags for AD, AC, and AN"
    bcftools +fill-tags "${OUTDIR}/datafiles/genotype_calls/${ID}_plinkfiltered".vcf.gz -Ob -o "{$ID}"_tagfilled.bcf -- -t AF,AC,AN
    bcftools index "{$ID}"_tagfilled.bcf;
done