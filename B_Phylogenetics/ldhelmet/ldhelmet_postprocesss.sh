#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script extracts estimates of recombination rates from binary ldhelmet output files.

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

printf "\n"
echo "|---------------Post-Processing ${CHR}---------------|"
printf "\n"

if [ -d "${RESULT_DIR}/${POP}" ]; then
        echo "${POP} result directory already exists."
    else
        echo "${POP} result directory does not exist. Creating it now..."
        mkdir -p "${RESULT_DIR}/${POP}"
    fi

# Have found that 20 threads on ocelote cluster can work for most max-likelihood analyses for bigger chromosomes. Configure your slurm script as needed

${PROGDIR}/LDhelmet/ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o ${RESULT_DIR}/${POP}/${POP}_${CHR}_STATS.txt ${INPUT_DIR}/${POP}/${CHR}.post

${PROGDIR}/LDhelmet/ldhelmet max_lk --num_threads 20 -l ${INPUT_DIR}/${POP}/${CHR}.lk -p ${INPUT_DIR}/${POP}/${CHR}.pade -s ${INPUT_DIR}/${POP}/${CHR}.ldhelmet.snps > ${RESULT_DIR}/${POP}/${POP}_${CHR}_maxlk.txt
