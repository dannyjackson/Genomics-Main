#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script generates demography estimates for a population using SMC++. This is a prerequisite for generating recombination maps using pyrho in A2.5.1_recombination_mapping.sh

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
Optional argument:
  -d  Boolean argument specifying whether to convert a copy of SMC++ output to demes YAML file (Default to false). Needed if want to use downstream Flexsweep analysis"
    exit 1
fi

DEMES=false

# Parse command-line arguments
while getopts ":p:d" option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        d) DEMES=true ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
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

if [ -d "${OUTDIR}/datafiles/demography/smc_inputs/${POP}" ];
        then
            echo "demography directory already exists, moving on!"
        else
            mkdir -p "${OUTDIR}/datafiles/demography/smc_inputs/${POP}"
fi

bcftools index ${VCF}

for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt);
do
    echo "Processing $chr into SMC++ Input"
    echo $POP:$POPSAMPLES
    COMMAND="$VCF ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/${POP}_${chr}.smc.gz $chr $POP:$POPSAMPLES --mask ${OUTDIR}/datafiles/snpable/${REF_ACC}_revised_mask.${chr}.mask.bed.gz"
    apptainer run ${PROGDIR}/smcpp_latest.sif vcf2smc $COMMAND
done

echo "Making SMC Input List"
ls ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/*.smc.gz > ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/smc_input_lst.txt
mapfile -t input_lst < ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/input_lst.txt
smc_inputs=$(echo "${smc_input_lst[@]}")


echo "Making and Plotting Model"
apptainer run ${PROGDIR}/smcpp_latest.sif estimate ${MUT_RATE} ${smc_inputs} -o ${OUTDIR}/datafiles/demography
apptainer run ${PROGDIR}/smcpp_latest.sif plot ${OUTDIR}/datafiles/demography/${OUTFNAME} ${OUTDIR}/datafiles/demography/model.final.json -c # Ensure to output csv file with model info for recombination mapping later
echo "Raw SMC++ Model outputted to model.final.json. CSV outputted to ${OUTFNAME}.csv"

if [ "$DEMES" = true ]; then
    echo "Converting model CSV file to Demes YAML file"
    python3 ${SCRIPTDIR}/Genomics-Main/general_scripts/convert_to_demes.py -c ${OUTDIR}/datafiles/demography/${OUTFNAME}.csv -t "years" -d "popsizes from smc++" -g 25 -o ${OUTFNAME}
    echo "Demes YAML outputted to ${OUTFNAME}.yaml"
fi