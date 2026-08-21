#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -m <recombination_map>"
    echo ""
    echo "This script phases an Individual-Scaffold VCF file (generated from A2.4) using BEAGLE5.5. This is to be used for the MSMC pipeline"
    echo "It is best run as a Slurm array that calls this script for each individual. Assuming that recombination maps are split by scaffold to reduce per-job computation cost"
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the individual to analyze."
    echo "  -m  Recombination Map file prefix (recommended). Assumes that map files are in datafiles/recombination_map directory and named in the format: PREFIX_SCAFFOLD.map"
    exit 1
}

# Parse command-line arguments
while getopts ":p:i:m:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) IND=${OPTARG} ;;
        m) MAPPREFIX=${OPTARG} ;;
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
echo "Recombination map prefix provided: ${MAPPREFIX}"

SCAFFOLD_LST=${OUTDIR}/referencelists/SCAFFOLDS.txt
while read -r SCAFFOLD; do

    VCF_IN="${OUTDIR}/datafiles/vcf/${IND}.${SCAFFOLD}.vcf"
    VCF_OUT="${OUTDIR}/datafiles/vcf2/${IND}.${SCAFFOLD}.phased"
    
    if [[ ! -f "$VCF_IN" ]]; then
        echo "Warning: Input VCF file for scaffold $SCAFFOLD not found. Exiting script. Please check path/filenames."
        exit 1
    fi
    if [[ -f "$VCF_OUT.vcf.gz" ]]; then
        echo "Phased VCF already exists for scaffold $SCAFFOLD. Skipping."
        continue
    fi

    if [[ -n "$MAPPREFIX" ]]; then
        MAP=${OUTDIR}/datafiles/recombination_map/${MAPPREFIX}_${SCAFFOLD}.map
        echo "Using Recombination map: ${MAP}"
        map_flag="map=${MAP}"
        if [[ ! -f "$MAP" ]]; then
            echo "Warning: Recombination Map not found at path: ${MAP}. Exiting script. Please check path/filenames"
            exit 1
        fi
    else
        echo "No recombination map provided. Proceeding without it. Note that BEAGLE will use a default recombination rate of 1"
        map_flag=""
    fi


    echo "Phasing VCF for scaffold $SCAFFOLD..."
    beagle gt=${VCF_IN} out=${VCF_OUT} ${map_flag}
    bcftools index ${VCF_OUT}.vcf.gz

done < "$SCAFFOLD_LST"

echo "VCF phasing for ${IND} completed."
date
