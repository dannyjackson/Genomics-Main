#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script generates demography estimates for a population using SMC++. This is a prerequisite for generating recombination maps using pyrho in A2.5.1_linkage_mapping.sh

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -i Individual name (can easily be passed through a slurm array)."
    exit 1
fi

DEMES=false

# Parse command-line arguments
while getopts p:s:r:v: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		s) SCAFFOLD=${OPTARG};;
        r) RUNPARAMS=${OPTARG};;
        v) VCF=${OPTARG};;
        d) DEMES=true
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"
source "${RUNPARAMS}"

printf "\n\n\n\n"
date
echo "Current script: A2.5.1_smc_demography.sh"

if [ -f "${OUTDIR}/datafiles/demography/" ];
        then
            echo "demography directory already exists, moving on!"
        else
            mkdir -p "${OUTDIR}/datafiles/demography/"
fi

apptainer run ${PROGDIR}/smcpp_latest.sif vcf2smc ${VCF} ${OUT_SMC_FILE} ${POPSET} --mask ${OUTDIR}/datafiles/snpable/${prefix}_revised_mask.bed.gz
apptainer run ${PROGDIR}/smcpp_latest.sif estimate ${MUT_RATE} ${OUT_SMC_FILE} -o ${OUTDIR}/datafiles/demography/ 
apptainer run ${PROGDIR}/smcpp_latest.sif plot ${OUTDIR}/datafiles/demography/${OUTFNAME} ${JSONMODEL} -c # Ensure to output csv file with model info for linkage mapping later

# Convert popsize csv to demes
if [ "$DEMES" = true ]; then
    python3 ${SCRIPT_DIR}/Genomics-Main/general_scripts/convert_to_demes.py -c ${OUTDIR}/datafiles/demography/${OUTFNAME}.csv -t 'years' -d "popsizes from smc++" -g 25 -o ${OUTFNAME}