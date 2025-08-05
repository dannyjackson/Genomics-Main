#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates input files for MSMC."
    echo "Required Argument:"
    echo "  -p   Path to parameter file (example in GitHub repository as params.sh)"
    echo "  -m   File Name of your unique project msmc params file"
    exit 1
fi

# Parse command-line arguments
while getopts "pm" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
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
# Check available modules (useful for debugging environment)
module list

# Source MSMC params file
source "${SCRIPTDIR}/${MSMCPARAMS}"


### Variables:
IND=`cat $1`
POP=$2

for s in `cat ${OUTDIR}/SCAFFOLDS.txt`
        do SCAFFOLD=$s

        MSMC_INPUT=${MSMCDIR}/input/msmc_input.${POP}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_multiInd"
        echo "Individuals: ${IND}"
        echo "Population: $POP"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"

        for ind in $IND
                do INDMASK=`ls ${MSMCDIR}/mask/ind/ind_mask.${ind}.${SCAFFOLD}.bed.gz`
                echo "--mask=$INDMASK " >> ${MSMCDIR}/mask/ind/${POP}.mask_file.$SCAFFOLD
                INDVCF=`ls ${MSMCDIR}/vcf2/${ind}.${SCAFFOLD}.phased.vcf.gz`
                echo $INDVCF >> ${MSMCDIR}/vcf2/${POP}.vcf_file.${SCAFFOLD}
        done

### Generate MSMC input files:
        if [ $METHOD == samtools ]
                then
                MASK_GENOME=${MSMCDIR}/mask/genom/${prefix}_revised_${SCAFFOLD}_mask.${k}.50.bed.gz

                echo "MAPPABILITY MASK: ${MASK_GENOME}"
                echo "Creating MSMC input file WITH individual mask (samtools)"
                #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
                ${MSMCTOOLS}/generate_multihetsep.py `cat ${MSMCDIR}/mask/ind/${POP}.mask_file.${SCAFFOLD}` --mask=$MASK_GENOME `cat ${MSMCDIR}/vcf2/${POP}.vcf_file.${SCAFFOLD}` > ${MSMC_INPUT} # without repeat mask

        elif [ $METHOD == gatk ]
                then
                echo "Creating MSMC input file WITHOUT individual mask (gatk)"
                MASK_GENOME=`ls ${MSMCDIR}/mask/${prefix}_${SCAFFOLD}.mask.${k}.50.bed.gz`
                #msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
                ${MSMCTOOLS}/generate_multihetsep.py --mask=$MASK_GENOME `cat ${MSMCDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}` > $MSMC_INPUT # without repeat mask

                # NOTE THAT THIS WAS CHANGED 10 FEB 2024
                # AND HAS NOT YET BEEN TESTED TO MAKE SURE IT WORKS
                                
                echo "Creating individual mask. Note that your input VCF should include ALL sites (variant & invariant)."
                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz
                                
                VCF_OUT=${VCF}.parsed.vcf
                ${MSMCTOOLS}/vcfAllSiteParser.py $SCAFFOLD $MASK_INDIV $VCF_OUT
                                
                echo "Creating MSMC input file with new individual mask"
                                
                ${MSMCTOOLS}/generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # with new repeat mask
        fi

done

echo "Done with script."
date

####
