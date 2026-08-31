#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates input files for MSMC."
    echo "Required Argument:"
    echo "  -p   Path to msmc parameter file (must source params_base.sh)"
    exit 1
fi

# Parse command-line arguments
while getopts ":p:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Ensure parameter file is provided and exists
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
elif [ ! -f "${PARAMS}" ]; then
    echo "Error: Parameter file '${PARAMS}' not found." >&2
    exit 1
fi

# Source the parameter file
source "${PARAMS}"

for SCAFFOLD in `cat ${OUTDIR}/referencelists/SCAFFOLDS.txt`
        do
        # Ensure list files are empty
        echo "" > ${OUTDIR}/datafiles/msmc/mask/ind/${POPNAME}.mask_file.${SCAFFOLD}
        echo "" > ${OUTDIR}/datafiles/msmc/vcf/${POPNAME}.vcf_file.${SCAFFOLD}

        MSMC_INPUT=${OUTDIR}/datafiles/msmc/input/msmc_input.${POPNAME}.${SCAFFOLD}.txt
        MASK_GENOME=${OUTDIR}/datafiles/snpable/${REF_ACC}_revised_mask.${SCAFFOLD}.mask.bed.gz

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_multiInd"
        echo "Individuals: ${INDFILE}"
        echo "Population: ${POPNAME}"
        echo "Scaffold: ${SCAFFOLD}"
        echo "MSMC input file: ${MSMC_INPUT}"
        echo "Genome Mask: ${MASK_GENOME}"

        for ind in $(cat ${INDFILE})
                do INDMASK=`ls ${OUTDIR}/datafiles/msmc/mask/ind/${ind}.${SCAFFOLD}.bed.gz`
                echo "--mask=$INDMASK " >> ${OUTDIR}/datafiles/msmc/mask/ind/${POPNAME}.mask_file.${SCAFFOLD}
                INDVCF=`ls ${OUTDIR}/datafiles/msmc/vcf/${ind}.${SCAFFOLD}.vcf.gz`
                echo $INDVCF >> ${OUTDIR}/datafiles/msmc/vcf/${POPNAME}.vcf_file.${SCAFFOLD}
        done

        # Generate MSMC input files:
        echo "Creating MSMC input file WITH individual mask"
        #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
        ${PROGDIR}/msmc-tools/generate_multihetsep.py `cat ${OUTDIR}/datafiles/msmc/mask/ind/${POPNAME}.mask_file.${SCAFFOLD}` \
            --mask=$MASK_GENOME `cat ${OUTDIR}/datafiles/msmc/vcf/${POPNAME}.vcf_file.${SCAFFOLD}` > ${MSMC_INPUT} # without repeat mask

done

echo "Done with script."
date

####
