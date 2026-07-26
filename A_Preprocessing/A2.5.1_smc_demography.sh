#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -v <vcf_file> -s <sample_list_file> -i <population_prefix> [-d]

This script generates demography estimates for a population using SMC++. This is a prerequisite for generating recombination maps using pyrho in A2.5.1_recombination_mapping.sh

Required argument:
  -p  Path to a parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -v  Path to Population VCF file
  -s  Path to Population Sample List File
  -i  Population prefix
Optional argument:
  -d  Boolean argument specifying whether to convert a copy of SMC++ output to demes YAML file (Default to false). Needed if want to use downstream Flexsweep analysis"
    exit 1
fi

DEMES=false

# Parse command-line arguments
while getopts ":p:v:s:i:d" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        v) VCF=${OPTARG} ;;
        s) SAMPLES=${OPTARG} ;;
        i) POP=${OPTARG} ;;
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

#bcftools index ${VCF}
# Check for bed mask file
#if [ -f "${OUTDIR}/datafiles/mask/${POP}_smc_mask.bed.gz" ];
#        then
#            echo "bed mask file already exists, moving on!"
#        else
#            # Get chrom order of VCF to ensure order in scaffold lengths file matches
#            zgrep -v '^#' ${VCF} | awk '{print $1}' | uniq > ${OUTDIR}/datafiles/genotype_calls/${POP}_vcf_chrom_order.txt
#            awk 'FNR==NR {pos[$1]=NR; next} $1 in pos {print $0, pos[$1]}' vcf_chrom_order.txt ${OUTDIR}/referencelists/scaffold_lengths.txt \
#            | sort -k3,3n \
#            | awk '{print $1"\t"$2}' > ${OUTDIR}/referencelists/scaffold_lengths.vcf_sorted.txt
#            # Generate bed file
#            bedtools complement -i ${VCF}.sorted -g ${OUTDIR}/referencelists/scaffold_lengths.vcf_sorted.txt > ${OUTDIR}/datafiles/mask/${POP}_smc_mask.bed
#            bgzip ${OUTDIR}/datafiles/mask/${POP}_smc_mask.bed
#            tabix ${OUTDIR}/datafiles/mask/${POP}_smc_mask.bed.gz
#fi

POPSAMPLES="$(cat ${SAMPLES} | paste -sd ",")"

for chr in $(cat ${OUTDIR}/referencelists/SCAFFOLDS.txt);
do
    echo "Processing $chr into SMC++ Input"
    COMMAND="$VCF ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/${POP}_${chr}.smc.gz $chr $POP:$POPSAMPLES --mask ${OUTDIR}/datafiles/mask/${POP}_sorted_mask.bed.gz"
    apptainer run ${PROGDIR}/smcpp_latest.sif vcf2smc $COMMAND --ignore-missing # After VCF filtering, some samples may be thrown out. Passing --ignore-missing to proceed with smc without them and log the sample names
done

echo "Making SMC Input List"
ls ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/*.smc.gz > ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/smc_input_lst.txt
mapfile -t smc_input_lst < ${OUTDIR}/datafiles/demography/smc_inputs/${POP}/smc_input_lst.txt
smc_inputs=$(echo "${smc_input_lst[@]}")


echo "Estimating Model"
apptainer run ${PROGDIR}/smcpp_latest.sif estimate ${MUT_RATE} ${smc_inputs} -o ${OUTDIR}/datafiles/demography

echo "Plotting Model"
apptainer run ${PROGDIR}/smcpp_latest.sif plot ${OUTDIR}/datafiles/demography/${POP} ${OUTDIR}/datafiles/demography/model.final.json -c # Also output model info to csv for recombination mapping later
echo "Raw SMC++ Model outputted to model.final.json. CSV outputted to ${POP}.csv"

if [ "$DEMES" = true ]; then
    echo "Converting model CSV file to Demes YAML file"
    python3 ${SCRIPTDIR}/Genomics-Main/general_scripts/convert_to_demes.py -c ${OUTDIR}/datafiles/demography/${POP}.csv -t "years" -d "SMC++ ${POP} popsizes" -g 5 -o ${POP}
    echo "Demes YAML outputted to ${POP}.yaml"
fi