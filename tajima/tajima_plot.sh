#!/bin/sh

# Tajima's D Computation Script

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <parameter_file> [-w <window_size>] [-s <step_size>] [-c <chrom_file>]

This script computes Tajima's D within a population of genomes using a genotype likelihood framework implemented in ANGSD.

It requires SAF files as input, which can be generated using the scripts available at:
    https://github.com/dannyjackson/Genomics-Main

It computes:
    - Genome-wide average Tajima's D
    - Sliding window Tajima's D
    - Per-SNP Tajima's D

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  Path to the parameter file (e.g., params_dxy.sh)
OPTIONAL ARGUMENTS:
    -w  Window size for sliding window analysis 
    -s  Step size for sliding window analysis
    -c  Path to a chromosome mapping file (optional)
EOF
    exit 1
fi


# Parse command-line arguments
while getopts "p:w:s:c:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        s) STEP=${OPTARG} ;;
        *) echo "Error: Invalid option '-${OPTARG}'" >&2; exit 1 ;;
    esac
done

# Ensure the parameter file is provided and exists
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
elif [ ! -f "${PARAMS}" ]; then
    echo "Error: Parameter file '${PARAMS}' not found." >&2
    exit 1
fi

# Source the parameter file
source "${PARAMS}"

# Ensure required variables are set
if [ -z "${OUTDIR}" ] || [ -z "${POP}" ] || [ -z "${WIN}" ] || [ -z "${STEP}" ]; then
    echo "Error: OUTDIR, POP, WIN, or STEP is not set in the parameter file." >&2
    exit 1
fi

WIN_OUT="${OUTDIR}/analyses/thetas/${POP}/${WIN}/${POP}.theta.thetasWindow.pestPG"

# Run R script for plotting
echo "Generating Manhattan plot from ${WIN_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}"

echo "Script completed successfully!"
