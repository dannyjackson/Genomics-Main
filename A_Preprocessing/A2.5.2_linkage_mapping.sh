#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script generates linkage map for a population using pyrho. It requires csv outputs from SMC++ in the datafiles/demography directory
Requires file with specific run parameters and sourcing from params_preprocessing.sh.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
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
echo "Current script: A2.5.2_linkage_mapping.sh"

if [ -f "${OUTDIR}/datafiles/linkage_map/" ];
        then
            echo "linkage map directory already exists, moving on!"
        else
        mkdir -p "${OUTDIR}/datafiles/linkage_map/"
fi

if [ -f "${OUTDIR}/datafiles/linkage_map/${POP}.hdf" ];
        then
            echo "lookup table already exists, moving on!"
        else
            pyrho make_table -n ${NUM_HAPS} -N ${N} --mu ${MUT_RATE} --logfile ${OUTDIR}/datafiles/linkage_map/pyrho_table_log.txt --outfile ${OUTDIR}/datafiles/linkage_map/${POP}.hdf --smcpp_file ${SMCFILE} --decimate_rel_tol 0.1
fi

# Run this to get probably better estimates of hyperparameters prior to running optimize
#pyrho hyperparam -n ${NUM_HAPS} --mu ${MUT_RATE} --blockpenalty 50,100 \
#	--windowsize 25,50 --logfile ${OUTDIR}/datafiles/linkage_map/ --tablefile ${OUTDIR}/datafiles/linkage_map/${POP}.hdf \
#	--num_sims 3 --num_threads ${THREADS}\
#	--smcpp_file ${SMCFILE} --outfile ${OUTDIR}/datafiles/linkage_map/${POP}_hyperparam_results.txt


for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt); do
    pyrho optimize --tablefile ${OUTDIR}/datafiles/linkage_map/${POP}.hdf \
        --vcffile ${VCFDIR}/${POP}_${chr}.vcf.gz \
        --outfile ${OUTDIR}/datafiles/linkage_map/${POP}_${chr}.rmap \
        --blockpenalty 50 --windowsize 50 \
        --logfile ${OUTDIR}/datafiles/linkage_map/pyrho_optimize_${chr}_log.txt \
        --fast-missing;
done

pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize ${NUM_HAPS} \
	--tablefile ${OUTDIR}/datafiles/linkage_map/${POP}.hdf \
	--outfile ${OUTDIR}/datafiles/linkage_map/${POP}_r2.txt