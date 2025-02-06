#!/bin/bash

# This script filters a vcf s.

# Function to display usage
usage() {
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates a mask for a reference genome to determine SNPable regions."
    echo
    echo "Required Arguments:"
    echo "  -p  Path to the parameter file (e.g., params.sh from the GitHub repository)"
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

source ${PARAMS}


# Filter VCF based on quality
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf" ]
        then
            echo "qualitysort vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by quality greater than 100"
            bcftools view -i 'QUAL>100' "${OUTDIR}/datafiles/genotype_calls/${ID}_snps_multiallelic.vcf" > "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf"

fi

# Filter VCF based on depth, remove indels
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered.recode.vcf" ]
        then
            echo "qualitysort vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by depth (min 2, max 8), remove indels"
            vcftools --vcf "${OUTDIR}/datafiles/genotype_calls/${ID}_qualitysort.vcf" --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered"
fi



# Further filtering using PLINK
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_mind2.vcf" ]
        then
            echo "plink filtered vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "filtering VCF file by geno 0.2, minimum ind 0.2, maf 0.01, remove indels"
            plink --vcf "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered.recode.vcf" \
            --allow-extra-chr --snps-only 'just-acgt' \
            --geno 0.02 --mind 0.2 --maf 0.01 \
            --recode vcf-iid --out "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_mind2"
fi

