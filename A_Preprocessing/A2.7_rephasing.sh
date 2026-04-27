#!/bin/bash

# THIS SCRIPT IS UNTESTED!!!!!!!!!!!!!!

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-b <bcf_file>]"
    echo ""
    echo "This script rephases low confidence variants in a BCF file with SAPPHIRE. PUMA CLUSTER REQUIRED FOR THIS SCRIPT."
    echo "It is best run as a Slurm array that calls this script for each individual/population." Meant to be used after stat-based phasing in step A2.6
    echo "Prior to running, you should split your phased BCF by chromosome and pass each split vcf in a slurm array."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -c  Name of chromosome present in VCF file."
    echo "  -b  Path to the BCF file to rephase. This file should be split by chromosome. NOTE: Use the bcftools_fill_tags.sh script as an example of converting VCF to BCF AND using +filltags bcftools plugin to populate allele frequency info (Required if BEAGLE phasing)"

    exit 1
}

# Parse command-line arguments
while getopts ":p:c:b:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        c) CHR=${OPTARG} ;;
        b) BCF=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$POP" || -z "$VCF" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Load parameters
source "${PARAMS}"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

echo Inputted BCF: $BCF

# Perform sapphire polishing

REPHASE_DIR="${OUTDIR}/datafiles/rephased_bcf"

if [ -d "${REPHASE_DIR}/bcfs" ]; then
    echo "Directory ${REPHASE_DIR} exists."
else
    echo "Directory ${REPHASE_DIR} does not exist. Creating it now."
    mkdir -p ${REPHASE_DIR}/bcfs
fi

echo "Polishing phased VCF for $CHR with SAPPHIRE..."

# Key File Paths
BCF_ANNOTATED="${REPHASE_DIR}/${CHR}_phased_PP_annotated.bcf"
BCF_REPHASED="${REPHASE_DIR}/bcfs/${CHR}_rephased.bcf"

EXTRACTED_PP="${REPHASE_DIR}/${CHR}_phased_PP_extract.bin"
ANNOTATED_CSV="${REPHASE_DIR}/${CHR}_phased_PP_annotated_samples.csv"
PP_INFO="${REPHASE_DIR}/${CHR}_pp_info.tsv"
HEADER="${REPHASE_DIR}/header.txt"
CRAM_PATH="${OUTDIR}/datafiles/indelrealignment/" # Should be a directory containing all individual CRAM files

# Extract all variants with AF < 0.01 and replace heterozygous by "0.5" and homozygous by "."
bcftools filter -i 'INFO/AF < 0.01' $BCF -Ou | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' | sed -e 's/0|0/./g' -e 's/0|1/0.5/g' -e 's/1|0/0.5/g' -e 's/1|1/./g' > $PP_INFO
# bgzip and tabix annotate the file
bgzip $PP_INFO
tabix -s1 -b2 -e2 $PP_INFO.gz

if [ -f "${HEADER}" ];
    then
        echo "Header file already exists, moving on!"
    else
        echo '##FORMAT=<ID=PP,Number=1,Type=Float,Description="Phasing Probability (PP) field">' > ${HEADER}
fi

bcftools annotate -a ${PP_INFO}.gz -h ${HEADER} -c CHROM,POS,ID,REF,ALT,FORMAT/PP $BCF -Ob -o $BCF_ANNOTATED

# Sanity check the annotated VCF
#bcftools view "$BCF_ANNOTATED" | less

${PROGDIR}/sapphire/pp_extractor/pp_extract -f ${BCF_ANNOTATED} --pp-from-maf --maf-threshold 0.01 -o ${EXTRACTED_PP}
# Copy the file as SAPPHIRE will modify it, with this we will be able to compare
cp "$EXTRACTED_PP" "${EXTRACTED_PP}.original"

bcftools query --list-samples $BCF_ANNOTATED | \
awk -v c=$CRAM_PATH '{ print NR-1 "," $0 "," c $0 ".cram" }' > $ANNOTATED_CSV

${PROGDIR}/sapphire/phase_caller/phase_caller \
-f $BCF_ANNOTATED -S $ANNOTATED_CSV --cram-path-from-samples-file \
-b "$EXTRACTED_PP" -t $THREADS

#./sapphire/bin_tools/bin_diff -a "$EXTRACTED_PP_ORIGINAL" \
#-b "$EXTRACTED_PP" \
#-f $VCF_OUT \
#-S $ANNOTATED_CSV \
#--extra-info --more | less

${PROGDIR}/sapphire/pp_update/pp_update -f $BCF_ANNOTATED -b $EXTRACTED_PP -o $BCF_REPHASED



#Remove large intermediate files
rm "$BCF_ANNOTATED"

echo "BCF phasing completed."
date