#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script uses FlexSweep to estimate various selection statistics and predict the likelihood of a selective sweep in the population.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).

Required argument:
  -p  Path to flexsweep parameter file (must source from params_base.sh)."
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
echo "Current script: flexsweep_1_analysis.sh"

if [ ! -d "${OUTDIR}/analyses/flexsweep/${POPNAME}" ]; then
  echo "Directory for flexsweep output for ${POPNAME} does not exist. Creating it now..."
  mkdir -p "${OUTDIR}/analyses/flexsweep/${POPNAME}" # -p creates parent directories if they don't exist
else
  echo "Directory for flexsweep output for ${POPNAME} already exists."
fi

echo "Starting Simulations"
flexsweep simulator --sample_size ${NUM_HAPS} --demes ${DEMES} --output_folder ${OUTDIR}/analyses/flexsweep/${POPNAME}  --nthreads ${THREADS} --num_simulations ${SIMULATIONS}

echo "Estimating feature vectors from simulations"
flexsweep fvs-discoal --simulations_path ${OUTDIR}/analyses/flexsweep/${POPNAME}  --nthreads ${THREADS}

echo "Estmating feature vectors from vcfs"
# BCF or GZipped VCF files required. Tabix (samtools) required for indexing.
# can comment out recombination map flag if do not have one
flexsweep fvs-vcf --vcf_path ${VCFDIR} --recombination_map ${REC_MAP} --nthreads ${THREADS} --suffix ${POPNAME}

echo "Starting CNN"
flexsweep cnn  --train_data ${OUTDIR}/analyses/flexsweep/${POPNAME}/fvs.parquet --predict_data ${VCFDIR}/fvs_${POPNAME}.parquet --output_folder ${OUTDIR}/analyses/flexsweep/${POPNAME}

echo "Done"