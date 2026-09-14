#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script uses FlexSweep to estimate various selection statistics and predict the likelihood of a selective sweep in the population.

Required argument:
  -p  Path to flexsweep parameter file (must source from params_base.sh).
Optional Argument:
  -r  Specify whether to use Recombination map in VCF feature estimation. Options include 'true' or 'false'. Defaults to true."
    exit 1
fi

USE_RECMAP=true

# Parse command-line arguments
while getopts p:r: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) USE_RECMAP=${OPTARG};;
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
echo "Current script: flexsweep_1_analysis.sh"

if [ ! -d "${OUTDIR}/analyses/flexsweep/${POPNAME}" ]; then
  echo "Directory for flexsweep output for ${POPNAME} does not exist. Creating it now..."
  mkdir -p "${OUTDIR}/analyses/flexsweep/${POPNAME}" # -p creates parent directories if they don't exist
else
  echo "Directory for flexsweep output for ${POPNAME} already exists."
fi

if [[ -d "${OUTDIR}/analyses/flexsweep/${POPNAME}/neutral" && "${OUTDIR}/analyses/flexsweep/${POPNAME}/sweep" ]];
    echo "Neutral and Sweep simulation directories found. Assuming simulations are done and feature vectors estimated. Moving on..."
else
    echo "Starting Simulations"
    flexsweep simulator --sample_size ${NUM_HAPS} --demes ${DEMES} --output_folder ${OUTDIR}/analyses/flexsweep/${POPNAME}  --nthreads ${THREADS} --num_simulations ${SIMULATIONS}

    echo "Estimating feature vectors from simulations"
    flexsweep fvs-discoal --simulations_path ${OUTDIR}/analyses/flexsweep/${POPNAME}  --nthreads ${THREADS}

echo "Estmating feature vectors from vcfs"
if [ "$USE_RECMAP" = "true" ]; then
    rec_map_flag="--recombination_map ${REC_MAP}"
else
    rec_map_flag=""
    echo "No recombination map passed"
# BCF or GZipped VCF files required. Tabix (samtools) required for indexing.
flexsweep fvs-vcf --vcf_path ${VCFDIR} ${rec_map_flag} --nthreads ${THREADS} --suffix ${POPNAME}

echo "Starting CNN"
flexsweep cnn  --train_data ${OUTDIR}/analyses/flexsweep/${POPNAME}/fvs.parquet --predict_data ${VCFDIR}/fvs_${POPNAME}.parquet --output_folder ${OUTDIR}/analyses/flexsweep/${POPNAME}

echo "Done"