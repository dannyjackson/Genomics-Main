#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -c <chromosome> -b <bcf_file> -m <minor_allele_freq_threshold> -o <output_dir>"
    echo ""
    echo "This script rephases low confidence variants in a BCF file with SAPPHIRE using Minor Allele Frequency metrics. This is meant to be used with BEAGLE-phased inputs."
    echo "It is best run as a Slurm array that calls this script for each chromosome within a population." Meant to be used after stat-based phasing in step A2.6
    echo "PUMA CLUSTER REQUIRED FOR THIS SCRIPT."
    echo "Prior to running, you should:
            1) Convert your phased VCF to BCF. NOTE: You must have allele frequency (AF) INFO field populated (A2.2).
            2) Split your phased BCF by chromosome.
            3) Convert your individual BAM files into CRAM files. This script assumes that they are in a directory called datafiles/crams"
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -c  Name of chromosome present in BCF file."
    echo "  -b  Path to the BCF file to rephase. This file should contain sites from only a single chromosome. NOTE: You must have allele frequency (AF) INFO field populated (A2.2)."
    echo "  -m  Minor Allele Frequnency Threshold. Defaults to 0.01"
    echo "  -o  Path to preferred outdirectory."

    exit 1
}

MAF=0.01

# Parse command-line arguments
while getopts ":p:c:b:m:o:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        c) CHR=${OPTARG} ;;
        b) BCF=${OPTARG} ;;
        m) MAF=${OPTARG} ;;
        o) REPHASE_DIR=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Load parameters
source "${PARAMS}"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

echo Inputted BCF: $BCF

# Perform sapphire polishing

if [ -d "${REPHASE_DIR}/bcfs" ]; then
    echo "Directory ${REPHASE_DIR}/bcfs exists."
else
    echo "Directory ${REPHASE_DIR}/bcfs does not exist. Creating it now."
    mkdir -p ${REPHASE_DIR}/bcfs
fi

echo "Polishing phased BCF for $CHR with SAPPHIRE..."

# Key File Paths
BCF_ANNOTATED="${REPHASE_DIR}/${CHR}_phased_PP_annotated.bcf"
BCF_REPHASED="${REPHASE_DIR}/bcfs/${CHR}_rephased.bcf"

EXTRACTED_PP="${REPHASE_DIR}/${CHR}_phased_PP_extract.bin"
ANNOTATED_CSV="${REPHASE_DIR}/${CHR}_phased_PP_annotated_samples.csv"
PP_INFO="${REPHASE_DIR}/${CHR}_pp_info.tsv"
HEADER="${REPHASE_DIR}/header.txt"
CRAM_PATH="${OUTDIR}/datafiles/crams" # Should be a directory containing all individual CRAM files

# Extract all variants with AF < inputted MAF and replace heterozygous by "0.5" and homozygous by "."
filter_str="INFO/AF<${MAF}"
bcftools filter -i $filter_str $BCF -Ou | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' | sed -e 's/0|0/./g' -e 's/0|1/0.5/g' -e 's/1|0/0.5/g' -e 's/1|1/./g' > $PP_INFO
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


echo "Extracting PP INFO from Minor Allele Frequency"
${PROGDIR}/sapphire/pp_extractor/pp_extract -f ${BCF_ANNOTATED} -o ${EXTRACTED_PP} --pp-from-maf --maf-threshold ${MAF} # the --pp-from-maf and --maf-threshold are required to handle non-SHAPEIT-phased BCFs.
# Copy the file as SAPPHIRE will modify it, with this we will be able to compare
cp "$EXTRACTED_PP" "${EXTRACTED_PP}.original"

echo "Saving PP INFO to CSV"
bcftools query --list-samples $BCF_ANNOTATED | \
awk -v c=$CRAM_PATH '{ print NR-1 "," $0 "," c $0 ".cram" }' > $ANNOTATED_CSV


echo "Beginning Phase Calling"
${PROGDIR}/sapphire/phase_caller/phase_caller \
-f $BCF_ANNOTATED -S $ANNOTATED_CSV --cram-path-from-samples-file \
-b "$EXTRACTED_PP" -t $THREADS

#./sapphire/bin_tools/bin_diff -a "$EXTRACTED_PP_ORIGINAL" \
#-b "$EXTRACTED_PP" \
#-f $VCF_OUT \
#-S $ANNOTATED_CSV \
#--extra-info --more | less

echo "Updating BCF file with new phasing info"
${PROGDIR}/sapphire/pp_update/pp_update -f $BCF_ANNOTATED -b $EXTRACTED_PP -o $BCF_REPHASED --no-pp # Using --no-pp to not update the Phasing probability field in the BCF (because it doesn't exist in BEAGLE outputs)

bcftools index $BCF_REPHASED

#Remove large intermediate files
rm "$BCF_ANNOTATED"

echo "BCF rephasing completed."
date