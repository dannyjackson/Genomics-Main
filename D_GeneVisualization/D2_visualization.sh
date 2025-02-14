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
    echo "  -p   Path to base parameter file (see example params_base.sh in the repository)\n"
    echo "  -f   File name of gene list, excluding extension (made by D2_uniquegenes.sh script)\n"

    exit 1
fi

# Parse arguments
while getopts "p:f:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        f) FILE=${OPTARG} ;;
    esac
done

source ${PARAMS}

Rscript "${SCRIPTDIR}/Genomics-Main/D_GeneVisualization/D2_visualization.r" ${OUTDIR} ${SCRIPTDIR} ${FILE}

