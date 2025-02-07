#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script phases a BAM file using a VCF and mask generated through the SNPable pipeline."
    echo "It is best run as a Slurm array that calls this script for each individual."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the individual to analyze."
    echo "  -b  Directory containing BAM files."
    exit 1
}

# Parse command-line arguments
while getopts ":p:i:b:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) IND=${OPTARG} ;;
        b) BAMDIR=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$IND" || -z "$BAMDIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Ensure script directory is set
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

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

# Ensure genome reference is set
if [ -z "$REF" ]; then
    echo "Error: REF variable is not defined in params.sh"
    exit 1
fi

echo "Processing individual: $IND"

# Ensure BAM file exists
BAMFILE="${BAMDIR}/${IND}.realigned.bam"
if [[ ! -f "$BAMFILE" ]]; then
    echo "Error: BAM file for $IND not found in $BAMDIR. Exiting."
    exit 1
fi

# Ensure required directories exist
mkdir -p "${OUTDIR}/vcf2" "${OUTDIR}/stats"

# Ensure the scaffold list exists
SCAFFOLD_LIST="${OUTDIR}/referencelists/SCAFFOLDS.txt"
if [[ ! -f "$SCAFFOLD_LIST" ]]; then
    echo "Error: Scaffold list file not found: $SCAFFOLD_LIST"
    exit 1
fi

# Iterate through scaffolds
while read -r SCAFFOLD; do
    echo "Processing scaffold: $SCAFFOLD"
    
    VCF_IN="${OUTDIR}/vcf/${IND}.${SCAFFOLD}.vcf"
    VCF_OUT="${OUTDIR}/vcf2/${IND}.${SCAFFOLD}.${phasing}.vcf.gz"

    if [[ ! -f "$VCF_IN" ]]; then
        echo "Warning: Input VCF file for scaffold $SCAFFOLD not found. Skipping."
        continue
    fi

    if [[ -f "$VCF_OUT" ]]; then
        echo "Phased VCF already exists for scaffold $SCAFFOLD; skipping."
        continue
    fi

    echo "Phasing VCF for scaffold $SCAFFOLD..."
    
    # Remove leading spaces from VCF
    sed -i.bak 's/^ //g' "$VCF_IN"

    # Run whatshap phasing
    whatshap phase --reference "$REF" --ignore-read-groups \
                   -o "$VCF_OUT" "$VCF_IN" "$BAMFILE"
    
    if [[ $? -ne 0 ]]; then
        echo "Error: whatshap phase failed for scaffold $SCAFFOLD."
        exit 1
    fi

    # Run whatshap stats
    STATS_OUT="${OUTDIR}/stats/${IND}.${SCAFFOLD}.${prefix}.minDP10.${phasing}.stats.tsv"
    whatshap stats --tsv="$STATS_OUT" "$VCF_OUT"

    if [[ $? -ne 0 ]]; then
        echo "Error: whatshap stats failed for scaffold $SCAFFOLD."
        exit 1
    fi

done < "$SCAFFOLD_LIST"

echo "VCF phasing completed."
date
