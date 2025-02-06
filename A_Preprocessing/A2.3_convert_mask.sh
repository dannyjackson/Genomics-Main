#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -p <parameter_file>"
    echo "This script converts a mask for a reference genome to determine which regions are SNPable."
    echo
    echo "Required Arguments:"
    echo "  -p    Path to the parameter file (example provided in the GitHub repository as params.sh)"
    exit 1
}

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    usage
fi

# Parse arguments
while getopts ":p:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        *) echo "Error: Invalid option -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure a parameter file is provided
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    usage
fi

# Load parameters
if [ ! -f "${PARAMS}" ]; then
    echo "Error: Parameter file '${PARAMS}' not found." >&2
    exit 1
fi

source "${PARAMS}"

# Ensure required environment variables are set
for var in PROGDIR PROJHUB OUTDIR REF; do
    if [ -z "${!var}" ]; then
        echo "Error: ${var} is not set in the parameter file." >&2
        exit 1
    fi
done

SNPABLE_SCRIPT_PATH="${PROGDIR}/seqbility-20091110"  # Directory with SNPable scripts

echo "Converting mask format..."
date

# Convert mask format
INPUT_MASK="${OUTDIR}/datafiles/snpable/${prefix}_mask.150.50.fa"
OUTPUT_MASK="${OUTDIR}/datafiles/snpable/${prefix}_revised_mask.150.50.fa"

if [ ! -f "${INPUT_MASK}" ]; then
    echo "Error: Input mask file '${INPUT_MASK}' not found." >&2
    exit 1
fi

sed 's/>>/>/g' "${INPUT_MASK}" > "${OUTPUT_MASK}"

echo "Creating mappability mask..."
date

# Create mappability mask
if [ ! -x "$(command -v python2)" ]; then
    echo "Error: python2 is not installed or not in PATH." >&2
    exit 1
fi

if [ ! -f "${PROGDIR}/msmc-tools/makeMappabilityMask.py" ]; then
    echo "Error: makeMappabilityMask.py not found in '${PROGDIR}/msmc-tools/'." >&2
    exit 1
fi

python2 "${PROGDIR}/msmc-tools/makeMappabilityMask.py"

echo "Mappability mask conversion completed."
date

