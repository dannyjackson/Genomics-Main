#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

Use GONE2 to generate estimates of recent Ne from variant data. All input files should be placed into a single input directory and named with the same prefix (as required by GONE2).

Required argument:
  -p  Path to the gone parameter file (must reference base parameter file)."
    exit 1
fi

REC_RATE=NONE
GENO_DTYPE=0

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
echo "Current script: gone_2_analysis.sh"

echo "Pulling input files from directory: ${INPUTDIR}"
echo "Using genotype data type option: ${GENO_DTYPE}"


if [ "$GENO_DTYPE" == "2" ]; then
    echo "Assuming Input file is VCF"
    INPUTFILE=$(find ${INPUTDIR}/*.vcf)
else
    echo "Assuming Input file is PLINK .ped"
    INPUTFILE=$(find ${INPUTDIR}/*.ped)
fi


RECOMB_OPT=""
if [ "$REC_RATE" == "NONE" ]; then
    echo "Assuming properly named recombination map provided in input directory..."
    echo "Running GONE with ${NUMIND} individuals with recombination map..."
else
    echo "User provided constant recombination rate ${REC_RATE}. Proceeding without recombination map"
    echo "Running GONE with ${NUMIND} individuals with constant recombination rate of ${REC_RATE}..."
    RECOMB_OPT="-r $RECOMB_RATE"
fi

RESULT_DIR=${OUTDIR}/analyses/gone2/${OUTNAME}
if [ ! -d "$RESULT_DIR" ]; then
  echo "Directory for gone2 output for ${OUTNAME} does not exist. Creating it now..."
  mkdir -p "$RESULT_DIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone2 output for ${OUTNAME} already exists. Moving on..."
fi

${PROGDIR}/GONE2/gone2 ${INPUTFILE} -g $GENO_DTYPE -i $NUMIND -t ${THREADS} $RECOMB_OPT -s ${NUM_SNPS} -S ${SEED}

mv ${INPUTFILE}_GONE2_Ne ${RESULT_DIR}/${OUTNAME}_GONE2_Ne
mv ${INPUTFILE}_GONE2_d2 ${RESULT_DIR}/${OUTNAME}_GONE2_d2
mv ${INPUTFILE}_GONE2_STATS ${RESULT_DIR}/${OUTNAME}_GONE2_STATS

echo "Completed GONE Analysis for $OUTNAME"