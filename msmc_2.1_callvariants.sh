#!/bin/bash

# Ensure script directory and output directory are defined
if [ -z "$SCRIPTDIR" ]; then
    echo "Error: SCRIPTDIR is not defined. Please set this variable."
    exit 1
fi

# Load parameters
if [ ! -f "${SCRIPTDIR}/params.sh" ]; then
    echo "Error: params.sh not found in ${SCRIPTDIR}"
    exit 1
fi
source "${SCRIPTDIR}/params.sh"

# Create the output directory for genotype calls if not present
mkdir -p "${OUTDIR}/datafiles/genotype_calls"
# Move to the output directory for genotype calls
cd "${OUTDIR}/datafiles/genotype_calls" || {
    echo "Error: Cannot change to output directory. Exiting."
    exit 1
}

# Generate VCF file with SNP calling
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_snps_multiallelic.vcf" ]
        then
            echo "snps_multiallelic vcf file is present in genotype_calls directory, assuming it is already generated and moving on!"
        else
            echo "generating VCF file using bcftools"
            bcftools mpileup -Ou -f "$REF" -a FORMAT/AD,DP,INFO/AD,SP "$BAMDIR"/*.final.bam | \
                bcftools call -mv -V indels > "${OUTDIR}/datafiles/genotype_calls/${ID}_snps_multiallelic.vcf"
fi



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


# Compress and index the cleaned VCF
if [ -f "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_cleaned_zip.vcf.gz" ]
        then
            echo "filtered and cleaned and zipped vcf file is present in genotype_calls directory, assuming it is already generated and indexed and moving on!"
        else
            echo "compressing and indexing the final vcf"
            bgzip -c "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_cleaned.vcf" > "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_cleaned_zip.vcf"
            bcftools index "${OUTDIR}/datafiles/genotype_calls/${ID}_filtered_cleaned_zip.vcf.gz"
fi
