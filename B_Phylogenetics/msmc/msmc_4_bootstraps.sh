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
    echo "  -i   Individual or Population Name"
    exit 1
fi

# Parse command-line arguments
while getopts "pmi" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
        i) POP_OR_IND=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&3; exit 1 ;;
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


# Verify that output location for msmc_outputs exists
if [ -d ${OUTDIR}/bootstrap/outputs/${IND} ]; then
    echo ${IND} 'bootstrap out directory exists'
else
    mkdir ${OUTDIR}/bootstrap/outputs/${IND}
    mkdir ${OUTDIR}/bootstrap/outputs/${IND}/log_and_loop
    echo ${IND} 'bootstrap out directory created'
fi

#input for the bootstrapping
BS_INPUT=`for s in $(cat ${OUTDIR}/SCAFFOLDS.txt); do find ${OUTDIR}/single_indv_data/input_single/ -maxdepth 1 -name "msmc_input.*${IND}.${s}*.txt"; done`

echo 'BS_INPUT'
echo '======================='
echo $BS_INPUT
echo '======================='

#output from the bootstrapping 
BS_OUTPUT=${OUTDIR}/bootstrap/${IND}.bootstrap
echo $BS_OUTPUT

#echo "generating bootstraps for ${IND}"
${MSMCTOOLS}/multihetsep_bootstrap.py --out_dir_prefix $BS_OUTPUT --files $BS_INPUT

cd ${OUTDIR}/bootstrap
ls -d *${IND}.bootstrap_* > ${OUTDIR}/bs_file_lists/${IND}.bs_file_list.txt

echo "Finished generating bootstraps for ${IND}"
