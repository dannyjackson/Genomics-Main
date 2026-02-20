#!/bin/sh

# all parameters come from the msmc_param control file
# make edits there before using this script!

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates input files for MSMC."
    echo "Required Argument:"
    echo "  -p   Path to msmc parameter file"
    echo "  -b   Bootstrap being run"
    echo "  -i   Individual or Population Name"
    exit 1
fi

# Parse command-line arguments
while getopts ":p:b:i:" opt; do
    case "${opt}" in
        p) PARAMS=${OPTARG} ;;
        b) BOOT=${OPTARG} ;;
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

# Source the parameter file
source "${PARAMS}"
# Check available modules (useful for debugging environment)
module list



# Run MSMC on ONE of your bootstrapped MSMC Inputs for a given pop/individual
MSMC_BS=$(find ${OUTDIR}/datafiles/bootstraps/${POP_OR_IND}/$BOOT -maxdepth 2 -name "bootstrap_multihetsep*.txt")
echo $MSMC_BS

MSMC_OUTPUT=${OUTDIR}/analyses/msmc/bootstraps/${POP_OR_IND}/msmc_output.$BOOT
echo $MSMC_OUTPUT

# Verify that output location for msmc_outputs exists
if [ -d ${MSMC_OUTPUT} ]; then
    echo ${POP_OR_IND} 'bootstrap output directory exists'
else
    mkdir -p ${MSMC_OUTPUT}/log_and_loop
    echo ${POP_OR_IND} 'bootstrap output directory created'
fi

echo "INDEX: ${INDEX}"

echo "running msmc2 on bootstraps for ${BOOT}"

msmc2_Linux -t $THREADS -p $P_PAR -i $NUM_OPT -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_BS
    
echo "done with msmc bootstraps for ${BOOT}"