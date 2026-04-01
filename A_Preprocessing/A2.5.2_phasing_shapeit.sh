#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script phases a VCF file with SHAPEIT5 using a VCF split by scaffold generated previously through the SNPable pipeline. It can also polish the phased output with SAPPHIRE, which is highly recommended"
    echo "It is best run as a Slurm array that calls this script for each individual."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the population to analyze."
    echo "  -b  Directory containing BAM files."
    echo "  -s  Optional flag to enable SAPPHIRE polishing of phased VCFs (recommended)."
    exit 1
}

SAPPHIRE=false

# Parse command-line arguments
while getopts ":p:i:b:s" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        b) BAMDIR=${OPTARG} ;;
        s) SAPPHIRE=true ;; # False by default, recommended to set to true
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$IND" || -z "$BAMDIR" ]]; then
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

# Ensure genome reference is set
if [ -z "$REF" ]; then
    echo "Error: REF variable is not defined in params.sh"
    exit 1
fi

echo "Processing individual: $IND"


# Ensure the scaffold list exists
SCAFFOLD_LIST="${OUTDIR}/referencelists/SCAFFOLDS.txt"
if [[ ! -f "$SCAFFOLD_LIST" ]]; then
    echo "Error: Scaffold list file not found: $SCAFFOLD_LIST"
    exit 1
fi

# Iterate through scaffolds
while read -r SCAFFOLD; do
    echo "Processing scaffold: $SCAFFOLD"
    
    VCF_IN="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.vcf.gz"
    VCF_OUT="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.phased.vcf.gz"

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
    #sed -i.bak 's/^ //g' "$VCF_IN"

    # Run shapeit phasing
    shapeit4 phase_common --input "$VCF_IN" --region 1 --output "$VCF_OUT" --thread $THREADS
    
    if [[ $? -ne 0 ]]; then
        echo "Error: shapeit phase failed for scaffold $SCAFFOLD."
        exit 1
    fi

    # Perform sapphire polishing

    if [ "$SAPPHIRE" = true ]; then
        echo "Polishing phased VCF for scaffold $SCAFFOLD with SAPPHIRE..."
        sapphire phase --bam "${BAMDIR}/${IND}.realigned.bam" --vcf "$VCF_OUT" --out "${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.sapphire.vcf.gz" --threads $THREADS
        
        if [[ $? -ne 0 ]]; then
            echo "Error: SAPPHIRE polishing failed for scaffold $SCAFFOLD."
            exit 1
        fi
        
        #Remove the unpolished phased VCF
        rm "$VCF_OUT"
        
    fi
    # Ensure BAM file exists
    #BAMFILE="${BAMDIR}/${IND}.realigned.bam"
    #if [[ ! -f "$BAMFILE" ]]; then
    #    echo "Error: BAM file for $IND not found in $BAMDIR. Exiting."
    #    exit 1
    #fi



done < "$SCAFFOLD_LIST"

echo "VCF phasing completed."
date
