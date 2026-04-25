#!/bin/bash


# THIS SCRIPT IS UNTESTED!!!!!!!!!!!!!!

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-s]"
    echo ""
    echo "This script phases a VCF file with BEAGLE5.5 using a VCF file containing all individuals in a desired population."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Path to VCF to phase."
    echo "  -o  Output VCF file path."
    echo "  -m  Path to recombination map file."
    exit 1
}

WINDOW=40
OVERLAP=2

# Parse command-line arguments
while getopts ":p:i:o:m:w:e:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) VCF_IN=${OPTARG} ;;
        o) VCF_OUT=${OPTARG} ;;
        m) MAP=${OPTARG} ;;
        w) WINDOW=${OPTARG} ;;
        c) OVERLAP=${OPTARG} ;;
        e) MARKERFILE=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

beagle gt=${VCF_IN} out=${VCF_OUT} map=${MAP} window=${WINDOW} overlap=${OVERLAP} excludemarkers=${MARKERFILE}
bcftools index ${VCF_OUT}

echo "VCF phasing completed."
date
