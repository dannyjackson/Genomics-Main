#!/bin/bash


# THIS SCRIPT IS UNTESTED!!!!!!!!!!!!!!

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-s]"
    echo ""
    echo "This script phases a VCF file with BEAGLE5.5 using a VCF file containing all individuals in a desired population."
    echo "It is best run as a Slurm array that calls this script for each individual/population."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -v  Path to VCF to phase."
    echo "  -o  Output VCF file path."
    echo "  -m  Path to recombination map file."
    exit 1
}

# Parse command-line arguments
while getopts ":p:v:o:m:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        v) VCF_IN=${OPTARG} ;;
        o) VCF_OUT=${OPTARG} ;;
        m) MAP=${OPTARG} ;;
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

beagle gt=${VCF_IN} out=${VCF_OUT} map=${OUTDIR}/datafiles/linkage_map/${POP}_all.rmap
bcftools index ${VCF_OUT}

echo "VCF phasing completed."
date
