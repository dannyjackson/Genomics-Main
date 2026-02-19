#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script estimates fine-scale recombination rates on a chromosome using LDhelmet.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency
This script assumes that the input files are generated using the ldhelmet_input_generation.sh from vcftools

Required argument:
  -p  Path to the parameter file (e.g., params_ldhelmet.sh in the GitHub repository).
  -d  Name of population to analyze. Should match name of directory containing input files if following full pipeline
  -c  Chromosome Name to process"
    exit 1
fi

# Parse command-line arguments
while getopts p:d:c: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		d) POP=${OPTARG};;
        c) CHR=${OPTARG};;
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
echo "Current script: ldhelmet_run_mcmc.sh"

printf "\n"
echo "|---------------Generating Haplotype Config for ${CHR}---------------|"
printf "\n"
ldhelmet find_confs --num_threads ${THREADS} -w ${WINDOW_SIZE} -o ${INPUT_DIR}/${POP}/${CHR}.conf ${INPUT_DIR}/${POP}/${CHR}.ldhelmet.snps


printf "\n"
echo "|---------------Generating Likelihood Lookup Table for ${CHR}---------------|"
printf "\n"
ldhelmet table_gen --num_threads ${THREADS} -t ${MUT_RATE} -r ${REC_RATE_GRID} -c ${INPUT_DIR}/${POP}/${CHR}.conf -o ${INPUT_DIR}/${POP}/${CHR}.lk 

printf "\n"
echo "|---------------Generating Pade Coefficients for ${CHR}---------------|"
printf "\n"
ldhelmet pade --num_threads ${THREADS} -t ${MUT_RATE} -x ${PADE_COEF} --defect_threshold ${DEFECT} -c ${INPUT_DIR}/${POP}/${CHR}.conf -o ${INPUT_DIR}/${POP}/${CHR}.pade


printf "\n"
echo "|---------------Running MCMC for ${CHR}---------------|"
printf "\n"

# Only add flags for custom files if provided
OPT_FLAGS=""
if [ -n "${MUT_MATRIX}" ]; then
    OPT_FLAGS="${OPT_FLAGS} -m ${MUT_MATRIX}"
fi
if [ -n "${ANC_PRIOR}" ]; then
    OPT_FLAGS="${OPT_FLAGS} -a ${ANC_PRIOR}"
fi

ldhelmet rjmcmc --num_threads ${THREADS} -l ${INPUT_DIR}/${POP}/${CHR}.lk -p ${INPUT_DIR}/${POP}/${CHR}.pade -s ${INPUT_DIR}/${POP}/${CHR}.ldhelmet.snps \
                -b ${BLOCK_PENALTY} --burn_in ${BURN_IN} -n ${ITERATIONS} -o ${INPUT_DIR}/${POP}/${CHR}.post ${OPT_FLAGS}


printf "\n"
echo "|---------------Finished running MCMC for ${CHR}---------------|"
printf "\n"