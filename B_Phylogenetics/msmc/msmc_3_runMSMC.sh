#!/bin/sh


# all parameters come from the msmc_param control file
# make edits there before using this script!

source /.bashrc
micromamba activiate msmc_env

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates input files for MSMC."
    echo "Required Argument:"
    echo "  -p   Path to parameter file (example in GitHub repository as params.sh)"
    echo "  -m   File Name of your unique project msmc params file"
    echo "  -i   Individual or Population Name"
    echo "  -t   Either "data" or "bootstraps". Doesn't change out MSMC is run, just determines where to look for input files and such. 
                        (Wanted single script to run MSMC from. If developing this workflow from the start, it would be implemented cleaner, but this works for now during testing)"
    exit 1
fi

# Parse command-line arguments
while getopts "p:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
        i) POP_OR_IND=${OPTARG} ;;
        t) RUN_TYPE=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&5; exit 1 ;;
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

RUN_NAME=msmc_${POP_OR_IND}_${DATE}
echo ${RUN_NAME}

### Report settings/parameters:
date
echo "Script: run_MSMC.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: "
echo "Scaffolds: SCAFS_INPUT_${POP_OR_IND}"
echo "Iterations: $NUM_OPT"


if [ $RUN_TYPE == data ]; then

                if [ $NR_IND == 1 ]; then

                        find ${MSMCDIR}/input/msmc_input.${POP_OR_IND}.*.txt -size 0 -delete
                        ls ${MSMCDIR}/input/msmc_input.${POP_OR_IND}.*.txt | grep -v $sex_chr > ${MSMCDIR}/input/SCAFS_INPUT_${POP_OR_IND}
                else
                        for i in `cat ${POP_OR_IND}_IND.txt`
                do echo $i
                IND=$i

                for s in `cat SCAFFOLDS.txt`
                        do echo $s
                        ls ${MSMCDIR}/input/msmc_input.${IND}.${s}.txt >> ${MSMCDIR}/input/SCAFS_INPUT_${POP_OR_IND}
                done
                done
                fi


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

                        msmc2_Linux -t $THREADS -p $P_PAR -i $NUM_OPT -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

                        mv $MSMC_OUTPUT*loop.txt ${MSMCDIR}/output/log_and_loop/
                        mv $MSMC_OUTPUT*log ${MSMCDIR}/output/log_and_loop/
                else
                        echo "Running MSMC for $NR_IND individuals"
                        MSMC_INPUT=`cat ${MSMCDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
                        MSMC_OUTPUT=${MSMCDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}
                        
                        msmc2_Linux -t $THREADS -p $P_PAR -i $NUM_OPT -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_INPUT

                        mv $MSMC_OUTPUT*loop.txt ${MSMCDIR}/output/log_and_loop/
                        mv $MSMC_OUTPUT*log ${MSMCDIR}/output/log_and_loop/
                fi

else
                for x in `cat ${MSMCDIR}/bs_file_lists/${POP_OR_IND}.bs_file_list.txt`; do
                        BS_INPUT_DIR=${MSMCDIR}/bootstrap/$x
                        MSMC_BS=$(find ${MSMCDIR}/bootstrap/$x -maxdepth 2 -name "bootstrap_multihetsep*.txt")
                        echo $MSMC_BS
                        
                        MSMC_OUTPUT=${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}/msmc_output.$x
                        echo $MSMC_OUTPUT
                        
                        echo "running msmc2 on bootstraps for $x"
                        msmc2_Linux -t $THREADS i- $NUM_OPT -o $MSMC_OUTPUT $MSMC_BS
                        
                        mv $MSMC_OUTPUT*loop.txt ${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}/log_and_loop/
                        mv $MSMC_OUTPUT*log ${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}/log_and_loop/


fi

echo "done running msmc2 for ${POP_OR_IND}"
date
