#!/bin/sh

# Tajima's D Computation Script

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <parameter_file> [-w <window_size>] [-s <step_size>] [-c <chrom_file>]

This script plots the output of a z-transformed windowed analysis of various metrics (dxy, fst, Tajima's D, etc.)

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  Path to the parameter file (e.g., params_dxy.sh)
    -w  Window size for sliding window analysis 
    -m  Name of metric (dxy, fst, Tajima)
EOF
    exit 1
fi


# Parse command-line arguments
while getopts "p:w:c:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        m) METRIC=${OPTARG} ;;
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
if [ -z "${OUTDIR}" ] || [ -z "${POP}" ] || [ -z "${WIN}" ] || [ -z "${METRIC}" ]; then
    echo "Error: OUTDIR, POP, WIN, or METRIC is not set in the parameter file." >&2
    exit 1
fi


WIN_OUT="${OUTDIR}/analyses/${METRIC}/${POP}/${POP}.${METRIC}_"${WIN}.Ztransformed.csv"

# Run R script for plotting
echo "Generating Manhattan plot from ${WIN_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}"

echo "Script completed successfully!"
