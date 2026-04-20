#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script generates demography estimates for a population using SMC++. This is a prerequisite for generating recombination maps using pyrho in A2.5.1_linkage_mapping.sh
Requires file with specific run parameters and sourcing from params_preprocessing.sh.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
Optional argument:
  -d  Boolean argument specifying whether to convert a copy of SMC++ output to demes YAML file (Default to false)"
    exit 1
fi

DEMES=false

# Parse command-line arguments
while getopts p:s:v:d option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
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

printf "\n\n\n\n"
date
echo "Current script: A2.5.1_smc_demography.sh"

if [ -f "${OUTDIR}/datafiles/demography/smc_inputs" ];
        then
            echo "demography directory already exists, moving on!"
        else
            mkdir -p "${OUTDIR}/datafiles/demography/smc_inputs"
fi


for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt);
do
    echo "Processing $chr into SMC++ Input"
    apptainer run ${PROGDIR}/smcpp_latest.sif vcf2smc ${VCF} ${OUTDIR}/datafiles/demography/smc_inputs/${CHR}.smc.gz ${CHR} ${POPSET} --mask ${OUTDIR}/datafiles/snpable/${prefix}_revised_mask.bed.gz
done

echo "Making and Plotting Model"
apptainer run ${PROGDIR}/smcpp_latest.sif estimate ${MUT_RATE} $(ls ${OUTDIR}/datafiles/demography/smc_inputs) -o ${OUTDIR}/datafiles/demography
apptainer run ${PROGDIR}/smcpp_latest.sif plot ${OUTDIR}/datafiles/demography/${OUTFNAME} ${OUTDIR}/datafiles/demography/model.final.json -c # Ensure to output csv file with model info for linkage mapping later

if [ "$DEMES" = true ]; then
    echo "Converting model CSV file to Demes YAML file
    python3 ${SCRIPT_DIR}/Genomics-Main/general_scripts/convert_to_demes.py -c ${OUTDIR}/datafiles/demography/${OUTFNAME}.csv -t 'years' -d "popsizes from smc++" -g 25 -o ${OUTFNAME}