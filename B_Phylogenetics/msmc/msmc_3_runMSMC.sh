#!/bin/sh

# all parameters come from the msmc_param control file
# make edits there before using this script!

source ~/.bashrc
micromamba activiate msmc_env

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates input files for MSMC."
    echo "Required Argument:"
    echo "  -p   Path to project base parameter file (example in GitHub repository as params.sh)"
    echo "  -m   File Name of your unique project msmc params file"
    echo "  -i   Individual or Population Name"
    exit 1
fi

# Parse command-line arguments
while getopts ":p:m:i:" opt; do
    case "${opt}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
        i) POP_OR_IND=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}"
           exit 1 
           ;;
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

# Source/list needed param files and modules
source "${PARAMS}"
source "${SCRIPTDIR}/${MSMCPARAMS}"
module list
RUN_NAME=msmc_${POP_OR_IND}_${DATE}


if [ $NR_IND == 1 ]; then

	find ${MSMCDIR}/input/msmc_input.${POP_OR_IND}.*.txt -size 0 -delete
	ls ${MSMCDIR}/input/msmc_input.${POP_OR_IND}.*.txt | grep -v $sex_chr > ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}
else

    for s in `cat /xdisk/mcnew/finches/ljvossler/finches/SCAFFOLDS.txt`
        do echo $s
        ls ${MSMCDIR}/multi_indv_data/input/msmc_input.${POP_OR_IND}.${s}.txt >> ${MSMCDIR}/multi_indv_data/input/SCAFS_INPUT_${POP_OR_IND}
    done
fi


### Report settings/parameters:
date
echo "Script: msmc_3_runMSMC.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Haplotype Indices Used: ${INDEX}"
echo "Population or Individuals ID: $POP_OR_IND"
echo "MSMC Inputs Read from: SCAFS_INPUT_${POP_OR_IND}"
echo "Iterations: ${NUM_OPT}"


if [ $NR_IND == 1 ]
        then
        echo "Running MSMC for one individual"
        MSMC_INPUT=`cat ${MSMCDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${MSMCDIR}/output/msmc_output.${RUN_NAME}

        if [ -f "${MSMCDIR}/input/SCAFS_INPUT_${POP_OR_IND}" ]
                then
                        echo "MSMC_INPUTS: SCAFS_INPUT_${POP_OR_IND}_noLG9"
                        echo "MSMC_OUTPUT: $MSMC_OUTPUT"
                else
                        echo "MSMC_INPUT does not exist! Exiting now"
                        exit 1
        fi

else
        echo "Running MSMC for $NR_IND individuals"
        MSMC_INPUT=`cat ${MSMCDIR}/multi_indv_data/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${MSMCDIR}/multi_indv_data/output_2/msmc_output.${RUN_NAME}

fi

echo "Output File Location: ${MSMC_OUTPUT}"


# Run MSMC now that all necessary params are set
msmc2_Linux -t $THREADS -p $P_PAR -i $NUM_OPT -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_INPUT 

mv $MSMC_OUTPUT*loop.txt ${MSMCDIR}/multi_indv_data/output_2/log_and_loop/
mv $MSMC_OUTPUT*log ${MSMCDIR}/multi_indv_data/output_2/log_and_loop/


echo "done running msmc2 for ${POP_OR_IND}"
date
