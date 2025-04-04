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
while getopts "p:w:s:m:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        s) STEP=${OPTARG} ;;
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
if [ -z "${OUTDIR}" ] || [ -z "${POP}" ] || [ -z "${WIN}" ] || [ -z "${STEP}" ]; then
    echo "Error: OUTDIR, POP, WIN, or STEP is not set in the parameter file." >&2
    exit 1
fi

# Print script metadata
echo -e "\n\n-------------------------"
date
echo "Running script: $(basename "$0")"
echo "-------------------------"

# Generate site frequency spectrum (SFS) file if it doesn't exist
SFS_FILE="${OUTDIR}/datafiles/safs/${POP}.sfs"
if [ -f "${SFS_FILE}" ]; then
    echo "Skipping SFS generation; file already exists: ${SFS_FILE}"
else
    echo "Generating SFS file..."
    "${ANGSD}/misc/realSFS" -P 24 "${OUTDIR}/datafiles/safs/${POP}.saf.idx" > "${SFS_FILE}"
fi

# Calculate per-site thetas if they don't exist
THETA_OUT="${OUTDIR}/analyses/thetas/${POP}/${POP}.thetas.idx"
if [ -f "${THETA_OUT}" ]; then
    echo "Skipping theta calculation; file already exists: ${THETA_OUT}"
else
    echo "Calculating thetas per site..."
    "${ANGSD}/misc/realSFS" saf2theta "${OUTDIR}/datafiles/safs/${POP}.saf.idx" \
        -outname "${OUTDIR}/analyses/thetas/${POP}/${POP}" \
        -sfs "${SFS_FILE}"
fi

# Estimate genome-wide Tajima's D if the output exists
TAJIMA_OUT="${OUTDIR}/analyses/thetas/${POP}/${POP}.thetas.idx.pestPG"
if [ -f "${TAJIMA_OUT}" ]; then
    echo "Skipping genome-wide Tajima's D estimation; file already exists: ${TAJIMA_OUT}"
else
    echo "Estimating genome-wide Tajima's D..."
    "${ANGSD}/misc/thetaStat" do_stat "${THETA_OUT}"
fi

WINDOW_OUT="${OUTDIR}/analyses/thetas/${POP}/${WIN}/${POP}.theta.thetasWindow"
WIN_OUT="${OUTDIR}/analyses/thetas/${POP}/${WIN}/${POP}.theta.thetasWindow.pestPG"

# Always run sliding window Tajima's D
if [ -f "${WIN_OUT}" ]; then
    echo "Skipping sliding-window Tajima's D estimation; file already exists: ${WIN_OUT}"
else
    echo "Estimating sliding window Tajima's D..."
    "${ANGSD}/misc/thetaStat" do_stat "${THETA_OUT}" \
        -win "${WIN}" -step "${STEP}" \
        -outnames "${WINDOW_OUT}"
fi


# Process chromosome mapping if a file is provided
if [ -n "${CHR_FILE}" ] && [ -f "${CHR_FILE}" ]; then
    echo "Processing chromosome mapping file: ${CHR_FILE}..."
    while IFS=',' read -r first second; do
        echo "Replacing '${second}' with '${first}' in ${WIN_OUT}..."
        sed -i.bak "s/${second}/${first}/g" "${WIN_OUT}"
    done < "${CHR_FILE}"
    rm -f "${WIN_OUT}.bak"
else
    echo "No valid chromosome mapping file provided. Skipping CHROM processing."
fi

# z transform windowed data
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}"

Z_OUT="${OUTDIR}/analyses/tajima/${POP}/${POP}.Tajima_${WIN}.Ztransformed.csv"
sed -i 's/\"//g' $Z_OUT

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "Tajima" "${POP}"

echo "Script completed successfully!"
