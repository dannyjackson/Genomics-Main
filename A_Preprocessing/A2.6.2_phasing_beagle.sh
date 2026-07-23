#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -m <recombination_map> -b"
    echo ""
    echo "This script phases an Individual-Scaffold VCF file (generated from A2.4) using BEAGLE5.5. This is to be used for the MSMC pipeline"
    echo "It is best run as a Slurm array that calls this script for each individual."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the individual to analyze."
    echo "  -m  Path to optional recombination map file. (recommended)"
    echo "  -b  Boolean parameter to convert phased VCF to BCF. Recommended if want to rephase VCF later using A2.7. (default to false)"
    exit 1
}

BCF=false

# Parse command-line arguments
while getopts ":p:i:m:b" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) IND=${OPTARG} ;;
        m) MAP=${OPTARG} ;;
        b) BCF=true ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$IND" ]]; then
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

echo "Processing individual: $IND"


if [[ -n "$MAP" ]]; then
    echo "Recombination map provided."
    map_flag="map=${MAP}"
else
    echo "No recombination map provided. Proceeding without it. Note that BEAGLE will use a default recombination rate of 1"
    map_flag=""
fi


SCAFFOLD_LST=${OUTDIR}/referencelists/SCAFFOLDS.txt
while read -r SCAFFOLD; do

    VCF_IN="${OUTDIR}/datafiles/vcf/${IND}.${SCAFFOLD}.vcf"
    VCF_OUT="${OUTDIR}/datafiles/vcf2/${IND}.${SCAFFOLD}.phased"
    
    if [[ ! -f "$VCF_IN" ]]; then
        echo "Warning: Input VCF file for scaffold $SCAFFOLD not found. Skipping."
        continue
    fi

    if [[ -f "$VCF_OUT" ]]; then
        echo "Phased VCF already exists for scaffold $SCAFFOLD; skipping."
        continue
    fi

    echo "Phasing VCF for scaffold $SCAFFOLD..."
    beagle gt=${VCF_IN} out=${VCF_OUT} ${map_flag}
    bcftools index ${VCF_OUT}.vcf.gz

    if [ "$BCF" = true ]; then
        echo "Converting phased ${IND}-${SCAFFOLD} VCF to BCF"
        if [ -d "${OUTDIR}/datafiles/vcf2/bcfs" ]; then
            echo "Directory ${OUTDIR}/datafiles/vcf2/bcfs exists."
        else
            echo "Directory ${OUTDIR}/datafiles/vcf2/bcfs does not exist. Creating it now."
            mkdir -p ${OUTDIR}/datafiles/vcf2/bcfs
        fi
        bcftools view ${VCF_OUT} -O b -o ${OUTDIR}/datafiles/vcf2/bcfs/${IND}.${SCAFFOLD}.phased.bcf
        bcftools index ${OUTDIR}/datafiles/vcf2/bcfs/${IND}.${SCAFFOLD}.phased.bcf
    fi

done < "$SCAFFOLD_LST"

echo "VCF phasing for ${IND} completed."
date
