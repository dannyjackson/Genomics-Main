#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script uses FlexSweep to estimate various selection statistics and predict the likelihood of a selective sweep in the population.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).

Required argument:
  -p  Path to the population parameter file (must source from params_base.sh)."
    exit 1
fi

# Parse command-line arguments
while getopts p: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
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
echo "Current script: flexsweep_2_analysis.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] |; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi

if [ ! -d "${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}" ]; then
  echo "Directory for flexsweep output for ${POPNAME} does not exist. Creating it now..."
  mkdir -p "${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}" # -p creates parent directories if they don't exist
else
  echo "Directory for flexsweep output for ${POPNAME} already exists. WARNING: Existing files in this directory may be overwritten."
fi

flexsweep simulator --sample_size ${NUM_HAPS} --demes ${OUTDIR}/datafiles/demography/${POPNAME}.yaml --output_folder ${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}  --nthreads ${THREADS} --num_simulation ${NUM_SIMULATIONS}

flexsweep fvs-discoal --simulations_path ${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}  --nthreads ${THREADS}

flexsweep fvs-vcf --vcf_path ${OUTDIR}/datafiles/rephased_vcf/ --recombination_map ${REC_MAP} --nthreads ${THREADS} --pop ${POPNAME}

flexsweep cnn  --train_data ${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}/fvs.parquet --predict_data ${OUTDIR}/datafiles/rephased_vcf/fvs_${POPNAME}.parquet --output_folder ${OUTDIR}/analyses/flexsweep_outputs/${POPNAME}