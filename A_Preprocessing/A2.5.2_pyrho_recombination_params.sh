#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_dir> -s <path_to_smc_file> -i <population_prefix> -n <num_haplotypes>

This script generates pyrho lookup tables and outputs stats on optimal hyperparams to use for recombination mapping (A2.5.3) Requires csv outputs from SMC++

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -s  Path to SMC++ output file
  -i  Population prefix
  -n  Maximum haplotypes present (2x sample size)"
    exit 1
fi

# Parse command-line arguments
while getopts ":p:s:i:n:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        s) SMCFILE=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        n) NUM_HAPS=${OPTARG} ;;
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
echo "Current script: A2.5.2_pyrho_recombination_params.sh"

if [ -d "${OUTDIR}/datafiles/recombination_map/" ];
        then
            echo "recombination map directory already exists, moving on!"
        else
        mkdir -p "${OUTDIR}/datafiles/recombination_map/"
fi

if [ -f "${OUTDIR}/datafiles/recombination_map/${POP}.hdf" ];
        then 
            echo "lookup table already exists, moving on!"
        else
            pyrho make_table -n ${NUM_HAPS} --approx -N ${NUM_HAPS} --mu ${MUT_RATE} --logfile ${OUTDIR}/datafiles/recombination_map/pyrho_table_log.txt \
            --outfile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf --smcpp_file ${SMCFILE} --decimate_rel_tol 0.1 --numthreads ${THREADS}
fi

# Run this to get probably better estimates of hyperparameters prior to running optimize (A2.5.3)
pyrho hyperparam -n ${NUM_HAPS} --mu ${MUT_RATE} --blockpenalty 25,50,75,100 \
	--windowsize 25,50,75,100 --logfile ${OUTDIR}/datafiles/recombination_map/ --tablefile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf \
	--num_sims 5 --num_threads ${THREADS} \
	--smcpp_file ${SMCFILE} --outfile ${OUTDIR}/datafiles/recombination_map/${POP}_hyperparam_results.txt


echo "Done"