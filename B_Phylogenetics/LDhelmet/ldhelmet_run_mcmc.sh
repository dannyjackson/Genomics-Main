#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script estimates fine-scale reocombination rates on a chromosome using LDhelmet.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency

Required argument:
  -p  Path to the parameter file (e.g., params_ldhelmet.sh in the GitHub repository).
  -c  Sequence data filename. This can be in either a fasta or alternative .sites format (can easily be passed through a slurm array)."
    exit 1
fi

# Parse command-line arguments
while getopts p:c: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		c) INPUTFILE=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

CHR=${INPUTFILE%.*} # Extract chromosome name from input filename (removing extension)

# Load parameters
source "${PARAMS}"

if [ ! -d "${INPUT_DIR}/${CHR}" ]; then
    echo "${CHR} INPUT_DIR already exists."
else
    echo "${CHR} INPUT_DIR does not exist. Creating it now..."
    mkdir -p "${INPUT_DIR}/${CHR}"
fi

if [ ! -d "${RESULT_DIR}/${CHR}" ]; then
    echo "${CHR} RESULT_DIR already exists."
else
    echo "${CHR} RESULT_DIR does not exist. Creating it now..."
    mkdir -p "${RESULT_DIR}/${CHR}"
fi

printf "\n\n\n\n"
date
echo "Current script: run_ldhelmet.sh"

printf "\n"
echo "|---------------Generating Haplotype Config for ${CHR}---------------|"
printf "\n"
ldhelmet find_confs --num_threads ${THREADS} -w ${WINDOW_SIZE} -o ${INPUT_DIR}/${CHR}.conf ${INPUT_DIR}/${CHR}/${INPUTFILE}


printf "\n"
echo "|---------------Generating Likelihood Lookup Table for ${CHR}---------------|"
printf "\n"
ldhelmet table_gen --num_threads ${THREADS} -t ${MUT_RATE} -r ${REC_RATE_GRID} -c ${INPUT_DIR}/${CHR}.conf -o ${INPUT_DIR}/${CHR}.lk 

printf "\n"
echo "|---------------Generating Pade Coefficients for ${CHR}---------------|"
printf "\n"
ldhelmet pade --num_threads ${THREADS} -t ${MUT_RATE} -x ${PADE_COEF} --defect_threshold ${DEFECT} -c ${INPUT_DIR}/${CHR}.conf -o ${INPUT_DIR}/${CHR}.pade


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

ldhelmet rjmcmc --num_threads ${THREADS} -l ${INPUT_DIR}/${CHR}.lk -p ${INPUT_DIR}/${CHR}.pade -s ${INPUT_DIR}/${CHR}/${INPUTFILE} \
                -b ${BLOCK_PENALTY} --burn_in ${BURN_IN} -n ${ITERATIONS} -o ${RESULT_DIR}/${CHR}.post ${OPT_FLAGS}

printf "\n"
echo "|---------------Finished running MCMC for ${CHR}---------------|"
printf "\n"