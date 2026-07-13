#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script generates recombination map for a population using pyrho. It requires csv outputs from SMC++

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -c  Optional boolean parameter to label, merge, and convert pyrho rmap to CM units. Designed to output in format for flexsweep selection analysis."
    exit 1
fi

CM_POSTPROCESSING=false

# Parse command-line arguments
while getopts p:l option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        c) CM_POSTPROCESSING=true
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
echo "Current script: A2.5.2_recombination_mapping.sh"

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
            pyrho make_table -n ${NUM_HAPS} -N ${N} --mu ${MUT_RATE} --logfile ${OUTDIR}/datafiles/recombination_map/pyrho_table_log.txt \
            --outfile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf --smcpp_file ${SMCFILE} --decimate_rel_tol 0.1
fi

# Run this to get probably better estimates of hyperparameters prior to running optimize
#pyrho hyperparam -n ${NUM_HAPS} --mu ${MUT_RATE} --blockpenalty 50,100 \
#	--windowsize 25,50 --logfile ${OUTDIR}/datafiles/linkage_map/ --tablefile ${OUTDIR}/datafiles/linkage_map/${POP}.hdf \
#	--num_sims 3 --num_threads ${THREADS}\
#	--smcpp_file ${SMCFILE} --outfile ${OUTDIR}/datafiles/linkage_map/${POP}_hyperparam_results.txt


for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt); do
    echo "Optimizing $chr map"
    pyrho optimize --tablefile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf \
        --vcffile ${VCFDIR}/${POP}_${chr}.vcf.gz \
        --outfile ${OUTDIR}/datafiles/recombination_map/${POP}_${chr}.rmap \
        --blockpenalty 50 --windowsize 50 \
        --logfile ${OUTDIR}/datafiles/recombination_map/pyrho_optimize_${chr}_log.txt \
        --ploidy 2;
done

echo "Computing Recombination Map Stats"
pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize ${NUM_HAPS} \
	--tablefile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf \
	--outfile ${OUTDIR}/datafiles/recombination_map/${POP}_r2.txt


if [ "$CM_POSTPROCESSING" = true ]; then
    echo "Starting conversion to CM units"

    echo "Generating converted pyrho rmaps"
    python3 ${SCRIPTDIR}/Genomics-Main/general_scripts/pyrho_to_cm.py --scaffolds ${OUTDIR}/referencelists/SCAFFOLDS.txt --popname ${POP} --indir ${OUTDIR}/datafiles/recombination_maps

    echo "Merging converted pyrho rmaps"
    ls ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/*.rmap > ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/pyrho_converted_lst.txt
    mapfile -t rmap_lst < ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/pyrho_converted_lst.txt
    cat $(echo "${rmap_lst[@]}") > ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/${POP}_merged.rmap

fi


echo "Done"