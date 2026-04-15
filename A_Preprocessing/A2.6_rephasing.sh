#!/bin/bash

# THIS SCRIPT IS UNTESTED!!!!!!!!!!!!!!

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-v <vcf_file>]"
    echo ""
    echo "This script rephases low confidence variants in a VCF file with SAPPHIRE. PUMA CLUSTER REQUIRED FOR THIS SCRIPT."
    echo "It is best run as a Slurm array that calls this script for each individual/population." Meant to be used after phasing in ste A2.5
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the population to analyze."
    echo "  -v  Path to the VCF file to rephase."

    exit 1
}

# Parse command-line arguments
while getopts ":p:i:v:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        v) VCF=${OPTARG} ;;
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

echo Inputted VCF: $VCF

# Perform sapphire polishing

REPHASE_DIR="${OUTDIR}/datafiles/rephased_vcf"

if [ -d "${REPHASE_DIR}" ]; then
    echo "Directory ${REPHASE_DIR} exists."
else
    echo "Directory ${REPHASE_DIR} does not exist. Creating it now."
    mkdir -p ${REPHASE_DIR}
fi

echo "Polishing phased VCF for $POP with SAPPHIRE..."

# Key File Paths
VCF_ANNOTATED="${REPHASE_DIR}/${POP}_phased_PP_annotated.vcf"
VCF_REPHASED="${REPHASE_DIR}/${POP}_rephased.vcf"

EXTRACTED_PP="${REPHASE_DIR}/${POP}_phased_PP_extract.bin"
EXTRACTED_PP_ORIGINAL="${REPHASE_DIR}/${POP}_phased_PP_extract.bin.original"
ANNOTATED_CSV="${REPHASE_DIR}/${POP}_phased_PP_annotated_samples.csv"
PP_INFO="${REPHASE_DIR}/pp_info.tsv"
HEADER="${REPHASE_DIR}/header.txt"
CRAM_PATH="${OUTDIR}/datafiles/indelrealignment/${POP}.realigned.bam"

# Extract all variants with AF < 0.01 and replace heterozygous by "0.5" and homozygous by "."
bcftools filter -i 'INFO/AF<0.01' $VCF -Ou | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' | sed -e 's/0|0/./g' -e 's/0|1/0.5/g' -e 's/1|0/0.5/g' -e 's/1|1/./g' > $PP_INFO
# bgzip and tabix annotate the file
bgzip $PP_INFO
tabix -s1 -b2 -e2 $PP_INFO.gz

if [ -f "${HEADER}" ];
    then
        echo "Header file already exists, moving on!"
    else
        echo '##FORMAT=<ID=PP,Number=1,Type=Float,Description="SHAPEIT Phasing Probability (PP) field">' > ${HEADER}
fi

bcftools annotate -a ${PP_INFO}.gz -h ${HEADER} -c CHROM,POS,ID,REF,ALT,FORMAT/PP \
$VCF -Ob -o $VCF_ANNOTATED

# Sanity check the annotated VCF
bcftools view "$VCF_ANNOTATED" | less

${PROGDIR}/sapphire/pp_extractor/pp_extract -f ${VCF_ANNOTATED} --pp-from-maf --maf-threshold ${MAF} -o ${EXTRACTED_PP}
# Copy the file as SAPPHIRE will modify it, with this we will be able to compare
cp "$EXTRACTED_PP" "$EXTRACTED_PP_ORIGINAL"

bcftools query --list-samples $VCF_ANNOTATED | \
awk '{print NR-1 "," $0 ",${CRAM_PATH}"}' > \
$ANNOTATED_CSV

${PROGDIR}/sapphire/phase_caller/phase_caller \
-f $VCF_ANNOTATED \
-S $ANNOTATED_CSV \
--cram-path-from-samples-file \
-b "$EXTRACTED_PP" \
-t $(nproc)

#./sapphire/bin_tools/bin_diff -a "$EXTRACTED_PP_ORIGINAL" \
#-b "$EXTRACTED_PP" \
#-f $VCF_OUT \
#-S $ANNOTATED_CSV \
#--extra-info --more | less

${PROGDIR}/sapphire/pp_update/pp_update -f $VCF_ANNOTATED \
-b $EXTRACTED_PP \
-o $VCF_REPHASED



#Remove intermediate VCFs
rm "$VCF_ANNOTATED"

echo "VCF phasing completed."
date