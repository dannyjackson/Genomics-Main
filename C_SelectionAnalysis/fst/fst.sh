#!/bin/sh

# FST computation script using ANGSD
# Computes Fst between two genome groups using genotype likelihoods

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $(basename $0) -p <param_file> [-w <window_size>] [-s <step_size>]"
    echo "\nThis script computes Fst between two genome groups using ANGSD."
    echo "It requires SAF files as input and produces genome-wide Fst estimates,"
    echo "as well as sliding window Fst statistics.\n"
    echo "Required Argument:"
    echo "  -p   Path to parameter file (see example params_fst.sh in the repository)\n"
    echo "Optional Arguments:"
    echo "  -w   Window size for Fst scans (default: 10,000)"
    echo "  -s   Step size for Fst scans (default: 10,000)"
    exit 1
fi

# Parse arguments
while getopts "p:w:s:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        w) WIN=${OPTARG:-10000} ;;
        s) STEP=${OPTARG:-10000} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

# Ensure required parameter file is provided
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

# Set defaults if not provided
WIN=${WIN:-10000}
STEP=${STEP:-10000}

# Read chromosome file 
CHROM=$(cat "$CHR_FILE")


# Print script info
echo "\nRunning FST computation script"
date
echo "Window size: $WIN | Step size: $STEP"

# Check and generate 2D SFS prior
MLS_FILE="${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml"
if [ -f "$MLS_FILE" ]; then
    echo "MLS file already exists, skipping SFS prior computation."
else
    echo "Computing 2D SFS prior..."
    ${ANGSD}/misc/realSFS "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" > "$MLS_FILE"
fi

# Compute FST index
FST_INDEX="${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx"
if [ -f "$FST_INDEX" ]; then
    echo "FST index already exists, skipping computation."
else
    echo "Computing FST index..."
    ${ANGSD}/misc/realSFS fst index "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" \
        -sfs "$MLS_FILE" -fstout "${OUTDIR}/analyses/fst/${POP1}_${POP2}"
fi

# Compute global FST estimate
GLOBAL_FST_FILE="${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt"
if [ -f "$GLOBAL_FST_FILE" ]; then
    echo "Global FST estimate already exists, skipping computation."
else
    echo "Computing global FST estimate..."
    echo -e "FST.Unweight\tFST.Weight" > "$GLOBAL_FST_FILE"
    ${ANGSD}/misc/realSFS fst stats "$FST_INDEX" >> "$GLOBAL_FST_FILE"
fi



# Compute sliding window FST
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
if [ -f "$WIN_OUT" ]; then
    echo "Sliding window FST output already exists, skipping computation."
else
    echo "Computing sliding window FST..."
    ${ANGSD}/misc/realSFS fst stats2 "$FST_INDEX" -win "$WIN" -step "$STEP" > "$WIN_OUT"
fi


# Replace chromosome names if conversion file is provided
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed "s/$second/$first/g" "$WIN_OUT" >> "${WIN_OUT}.chrom.txt" 
    done < "$CHROM"
fi

# replace header (for whatever reason, it lacks a label for the fst column)
sed -i '1s/^region\tchr\tmidPos\tNsites$/region\tchr\tmidPos\tNsites\tfst/' "${WIN_OUT}.chrom.txt"


# z transform windowed data
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"
