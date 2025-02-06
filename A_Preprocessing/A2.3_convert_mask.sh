#!/bin/bash

# MSMC Pipeline Script
if [ $# -lt 1 ]
  then
    echo " This script generates a mask for a reference genome to determine which regions are "snpable".

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as params.sh"

  else
    while getopts p: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;

    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"

# Ensure required environment variables are set
if [ -z "$PROGDIR" ] || [ -z "$PROJHUB" ] || [ -z "$OUTDIR" ]; then
    echo "Error: PROGDIR, PROJHUB, or OUTDIR is not set." >&2
    exit 1
fi

snpable_script_path="${PROGDIR}/seqbility-20091110" # Directory with snpable scripts

# Check for required arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 -p <parameter_file> "
    exit 1
fi

while getopts ":p:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Ensure required arguments are provided
if [ -z "${PARAMS}" ] " ]; then
    echo "Error: Missing required arguments." >&2
    exit 1
fi

# Load parameters from external file
source "${PARAMS}"

echo "Converting mask format"
date

# Convert mask format
sed 's/>>/>/g' "${REF}_mask.150.50.fa" > "${REF}_revised_mask.150.50.fa"

echo "Creating mappability mask"
# Create mappability mask
python2 "${PROGDIR}/msmc-tools/makeMappabilityMask.py"

echo "Mappability mask conversion completed."
