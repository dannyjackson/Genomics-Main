#!/bin/sh

# Signasel script to detect significant allele frequency changes

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $(basename $0) -s <species> -o <output_directory> [-s <step_size>]"
    echo "\nThis script identifies sites that show significant allele frequency changes due to selection using Signasel. It takes as input two mafs files output by ANGSD."
    echo "Required Argument:"
    echo "  -s   Name of species"
    echo "  -o   Path to directory containing gzipped mafs files in format <species>_pre.mafs.gz and <species>_post.mafs.gz. This will also serve as the output directory."
    exit 1
fi

# Parse arguments
while getopts "s:o:" option; do
    case "${option}" in
        s) SP=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

module load R

Rscript ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/deltaAF/signasel_wrapper.r ${SP} ${OUTDIR}
