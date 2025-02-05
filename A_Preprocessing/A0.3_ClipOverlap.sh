#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script clips overlapping read pairs using BamUtils

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -i Individual name (can easily be passed through a slurm array)."
    exit 1
fi

# Parse command-line arguments
while getopts p:i: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		i) IND=${OPTARG};;
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
echo "Current script: A2_ClipOverlap.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$FASTAS" ] || [ -z "$BAMUTILBAM" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


# Clip overlapping read pairs using bamutils
echo "clipping" "${IND}" >> ${OUTDIR}/datafiles/clipoverlap/clippingstats.txt 

${BAMUTILBAM} clipOverlap --in ${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam --out ${OUTDIR}/datafiles/clipoverlap/${IND}.all.sorted.marked.clipped.bam --stats --params

echo "done " ${IND} >>${OUTDIR}/datafiles/clipoverlap/clippingstats.txt 
