#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_file> -s <sample_list_file> -i <population_prefix> -d <distinguished_lineages> -m <mask_file> [-c]

This script generates demography estimates for a population using SMC++. This is a prerequisite for generating recombination maps using pyrho in A2.5.1_recombination_mapping.sh
The model csv output can be converted to demes format for other analyses/preprocessing using convert_to_demes.py in general_scripts.

Required argument:
  -p  Path to a parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -v  Path to Population VCF file (This script assumes you are passing a full (non-chromosome-split) VCF and will loop through each chromosome separately.)
  -s  Path to Population Sample List File
  -i  Population prefix
Optional argument:
  -d  Specifying pair of distinguished lineages. If left blank, SMC++ will default to first samples in list. Should be a space separated string of two sample names
  -m  Path to a bgzipped and indexed bed file containing uncalled/lowcoverage regions not present in the VCF file (recommended)
  -c Optional boolean flag to do cv (SMC cross-validation model estimation)"
    exit 1
fi

DIST=""
MASK=""
CV=false

# Parse command-line arguments
while getopts ":p:v:s:i:d:m:c" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        v) VCF=${OPTARG} ;;
        s) SAMPLES=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        d) DIST=${OPTARG} ;;
        m) MASK=${OPTARG} ;;
        c) CV=true ;;
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

if [ -d "${OUTDIR}/datafiles/demography/${POP}/smc_inputs/" ];
        then
            echo "demography directory for ${POP} already exists, moving on!"
        else
            mkdir -p "${OUTDIR}/datafiles/demography/${POP}/smc_inputs/"
fi

bcftools index ${VCF}
# Check for bed mask file
if [[ -n "$MASK" ]]; then
    echo "Mask File provided in path: ${MASK}."
    mask_flag="--mask ${MASK}"
else
    echo "No mask for uncalled regions provided. Proceeding without it, but note that you probably want a mask to account for uncalled regions."
    mask_flag=""
fi

# Check for dstinguished lineages
if [[ -n "$DIST" ]]; then
    echo "Distinguished Lineages provided as follows: ${DIST}"
    dist_flag="-d ${DIST}"
else
    echo "No distinguished lineages provided. Proceeding without it. SMC will choose the first samples in your sample list."
    dist_flag=""
fi

POPSAMPLES="$(cat ${SAMPLES} | paste -sd ",")"

cd ${OUTDIR}/datafiles/demography/${POP} # Move to output directory here so that iterate.dat intermediate file is stored here.

for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt);
do
    if [ -f "${OUTDIR}/datafiles/demography/${POP}/smc_inputs/${POP}_${chr}.smc.gz" ]
        then
            echo "SMC input file for ${POP} ${chr} found, skipping..."
        else
            echo "Processing $chr into SMC++ Input"
            COMMAND="$VCF ${OUTDIR}/datafiles/demography/${POP}/smc_inputs/${POP}_${chr}.smc.gz $chr $POP:$POPSAMPLES ${dist_flag}" #${mask_flag}"
            apptainer run ${PROGDIR}/smcpp_latest.sif vcf2smc $COMMAND -c 50000 --ignore-missing # After VCF filtering, some samples may be thrown out. Passing --ignore-missing to proceed with smc without them and log the sample names
    fi
done

echo "Making SMC Input List"
ls ${OUTDIR}/datafiles/demography/${POP}/smc_inputs/*.smc.gz > ${OUTDIR}/datafiles/demography/${POP}/smc_inputs/smc_input_lst.txt
mapfile -t smc_input_lst < ${OUTDIR}/datafiles/demography/${POP}/smc_inputs/smc_input_lst.txt
smc_inputs=$(echo "${smc_input_lst[@]}")

if [ "$CV" = true ]
    then
        echo "Estimating Model using cross validation strategy"
        apptainer run ${PROGDIR}/smcpp_latest.sif cv ${MUT_RATE} ${smc_inputs} -o ${OUTDIR}/datafiles/demography/${POP} --folds 5 --cores ${THREADS} --thinning 3000
        SMC_OUTDIR=${OUTDIR}/datafiles/demography/${POP}/fold0
    else
        echo "Estimating Model using standard estimation strategy"
        apptainer run ${PROGDIR}/smcpp_latest.sif estimate ${MUT_RATE} ${smc_inputs} -o ${OUTDIR}/datafiles/demography/${POP} --cores ${THREADS} --thinning 3000
        SMC_OUTDIR=${OUTDIR}/datafiles/demography/${POP}
fi

echo "Plotting Model"
apptainer run ${PROGDIR}/smcpp_latest.sif plot ${SMC_OUTDIR} ${SMC_OUTDIR}/model.final.json -c # Also output model info to csv for recombination mapping later
echo "Raw SMC++ Model outputted to model.final.json. CSV outputted to ${POP}.csv"

echo "Done"