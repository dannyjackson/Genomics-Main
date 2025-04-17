#!/bin/sh

# FST plotting and outlier identification script for the output of depth and mapping filtered windowed FST analyses


# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $(basename $0) -p <param_file> [-w <window_size>] [-s <step_size>]"
    echo "\nThis script computes Fst between two genome groups using ANGSD."
    echo "It requires SAF files as input and produces genome-wide Fst estimates,"
    echo "as well as sliding window Fst statistics.\n"
    echo "Required Argument:"
    echo "  -p   Path to parameter file (see example params_fst.sh in the repository)\n"
    echo "  -f   Path to filtered FST file"
    exit 1
fi

# Parse arguments
while getopts "p:f:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        f) FILE=${OPTARG} ;;
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

# Read chromosome file 
CHROM=$(cat "$CHR_FILE")


# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    cp ${FILE} "${FILE}.numchrom" 
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "${FILE}.numchrom" 
    done < "$CHR_FILE"
fi



# z transform windowed data
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.filteredfiles.r" \
    "${OUTDIR}" "${CUTOFF}" "${FILE}.numchrom" 

Z_OUT="${FILE}.numchrom.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.filteredfiles.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "fst" 

echo "Script completed successfully!"
