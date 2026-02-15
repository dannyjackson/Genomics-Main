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
    echo "  -p   Path to parameter file (example in GitHub repository as params.sh)"
    echo "  -m   File Name of your unique project msmc params file"
    echo "  -b   Bootstrap being run"
    echo "  -i   Individual or Population Name"
    exit 1
fi

# Parse command-line arguments
while getopts ":p:m:b:i:" opt; do
    case "${opt}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
        b) BOOT=${OPTARG} ;;
        i) POP_OR_IND=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}"
           exit 1 
           ;;
    esac
done

echo $PARAMS
echo $MSMCPARAMS
echo $POP_OR_IND
echo $BOOT

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



# Run MSMC on ONE of your bootstrapped MSMC Inputs for a given pop/individual
MSMC_BS=$(find ${MSMCDIR}/bootstrap/$BOOT -maxdepth 2 -name "bootstrap_multihetsep*.txt")
echo $MSMC_BS

MSMC_OUTPUT=${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}/msmc_output.$BOOT
echo $MSMC_OUTPUT

echo "INDEX: ${INDEX}"

echo "running msmc2 on bootstraps for ${BOOT}"

msmc2_Linux -t $THREADS -p $P_PAR -i $NUM_OPT -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_BS
    
echo "done with msmc bootstraps for ${BOOT}"