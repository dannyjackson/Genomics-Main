#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_dir> -i <population_prefix> -b <block_penalty> -w <windowsize> [-c]

This script generates recombination map for a population using pyrho.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -v  Path to chromosome-split vcfs
  -i  Population prefix
  -b  Block Penalty size (Defaults to 50)
  -w  Window Size (Defaults to 50)
Optional argument:
  -c  Boolean parameter to label and convert pyrho rmap to CM units. Designed to output in format for flexsweep selection analysis."
    exit 1
fi

BLOCK=50
WINDOW=50
CM_POSTPROCESSING=false

# Parse command-line arguments
while getopts ":p:v:i:b:w:c" option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCFDIR=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        b) BLOCK=${OPTARG} ;;
        w) WINDOW=${OPTARG} ;;
        c) CM_POSTPROCESSING=true ;;
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
echo "Current script: A2.5.3_pyrho_recombination_mapping.sh"

for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt); do
    echo "Optimizing $chr map"
    pyrho optimize --tablefile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf \
        --vcffile ${VCFDIR}/${POP}_${chr}.vcf.gz \
        --outfile ${OUTDIR}/datafiles/recombination_map/${POP}_${chr}.rmap \
        --blockpenalty ${BLOCK} --windowsize ${WINDOW} \
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
    python3 ${SCRIPTDIR}/Genomics-Main/general_scripts/pyrho_to_cm.py --scaffolds ${OUTDIR}/referencelists/SCAFFOLDS.txt --pop ${POP} --indir ${OUTDIR}/datafiles/recombination_maps

    echo "Merging converted pyrho rmaps"
    ls ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/*.rmap > ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/pyrho_converted_lst.txt
    mapfile -t rmap_lst < ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/pyrho_converted_lst.txt
    cat $(echo "${rmap_lst[@]}") > ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/${POP}_merged.rmap
    echo Merged pyrho rmap saved to: ${OUTDIR}/datafiles/recombination_map/pyrho_cm_converted/${POP}_merged.rmap

fi


echo "Done"