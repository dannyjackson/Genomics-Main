#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_dir> -i <population_prefix> -b <block_penalty> -w <windowsize> [-c]

This script generates recombination map for a population using pyrho.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -v  Path to chromosome-split vcfs
  -i  Population prefix
  -b  Block Penalty size (Defaults to 50)
  -w  Window Size (Defaults to 50)"
    exit 1
fi

BLOCK=50
WINDOW=50

# Parse command-line arguments
while getopts ":p:v:i:b:w:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        v) VCFDIR=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        b) BLOCK=${OPTARG} ;;
        w) WINDOW=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"
module load parallel

printf "\n\n\n\n"
date
echo "Current script: A2.5.3_pyrho_recombination_mapping.sh"

CHROMS=${OUTDIR}/referencelists/SCAFFOLDS.txt


echo "Parallel map optimization..."
echo "Processing ${THREADS} chromosomes at a time..."
cat ${CHROMS} | parallel --jobs "$THREADS" "
    pyrho optimize --tablefile ${OUTDIR}/datafiles/recombination_map/${POP}.hdf \
    --vcffile ${VCFDIR}/${POP}_{}.vcf.gz \
    --outfile ${OUTDIR}/datafiles/recombination_map/${POP}_{}.rmap \
    --blockpenalty ${BLOCK} --windowsize ${WINDOW} \
    --logfile ${OUTDIR}/datafiles/recombination_map/${POP}_pyrho_optimize_{}_log.txt \
    --ploidy 2"

echo "Done"