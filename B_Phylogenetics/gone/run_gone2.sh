#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This uses GONE2 to generate estimates of recent Ne from plink (.ped/.map)files.
I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency

Required argument:
  -p  Path to the GONE parameter file (e.g., params_gone2.sh in the GitHub repository).
  -s  Population specific-parameter file.
  -r  Run name, required for providing a unique name to output files (especially when using in an array)."
    exit 1
fi

# Parse command-line arguments
while getopts p:s:r: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        s) POPPARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"
source "${POPPARAMS}"

printf "\n\n\n\n"
date
echo "Current script: run_gone2.sh"


echo "Running GONE for $POPNAME with $NUMIND individuals with recombination rate of $RECOMB_RATE..."
echo "NOTE: We are using Input files from: $INDIR"

RESULT_DIR=${OUTDIR}/analyses/gone2_outputs/${RUNNAME}/${POPNAME}

if [ ! -d "$RESULT_DIR" ]; then
  echo "Directory for gone2 output for ${RUNNAME} with ${POPNAME} does not exist. Creating it now..."
  mkdir -p "$RESULT_DIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone2 output for ${RUNNAME} with ${POPNAME} already exists. WARNING: Existing files in this directory may be overwritten."
fi

${PROGDIR}/GONE2/gone2 $INDIR/$POPNAME.ped -g $GENO_DTYPE -r $RECOMB_RATE -i $NUMIND -t 4 -o $RESULT_DIR/$POPNAME $EXTRA_FLAGS



echo "Completed GONE Analysis for $POPNAME"