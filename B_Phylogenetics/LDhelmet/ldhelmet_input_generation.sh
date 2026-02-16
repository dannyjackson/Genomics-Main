#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This uses VCFtools to generate input files for LDHelmet

Required argument:
  -p  Path to the parameter file (e.g., params_ldhelmet.sh in the GitHub repository).
  -v  VCF file path."
    exit 1
fi

# Parse command-line arguments
while getopts p:v: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCF=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

printf "\n\n\n\n"
date
echo "Current script: ldhelmet_input_generation.sh"

if [ ! -d "${INPUT_DIR}/${CHR}" ]; then
    echo "${CHR} INPUT_DIR already exists."
else
    echo "${CHR} INPUT_DIR does not exist. Creating it now..."
    mkdir -p "${INPUT_DIR}/${CHR}"
fi