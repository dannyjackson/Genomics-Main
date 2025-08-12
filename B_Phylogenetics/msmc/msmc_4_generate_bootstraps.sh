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

echo $PARAMS
echo $MSMCPARAMS
echo $POP_OR_IND

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
echo MSMCDIR: ${MSMCDIR}
echo OUTDIR: ${OUTDIR}


# Verify that output location for msmc_outputs exists
if [ -d ${MSMCDIR}/bootstrap/outputs/${POP_OR_IND} ]; then
    echo ${IND} 'bootstrap out directory exists'
else
    mkdir ${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}
    mkdir ${MSMCDIR}/bootstrap/outputs/${POP_OR_IND}/log_and_loop
    echo ${POP_OR_IND} 'bootstrap out directory created'
fi

#input for the bootstrapping
BS_INPUT=`for s in $(cat ${OUTDIR}/SCAFFOLDS.txt); do find ${MSMCDIR}/multi_indv_data/input/ -maxdepth 1 -name "msmc_input.*${POP_OR_IND}.${s}*.txt"; done`

echo 'BS_INPUT'
echo '======================='
echo $BS_INPUT
echo '======================='

#output from the bootstrapping 
BS_OUTPUT=${MSMCDIR}/bootstrap/${POP_OR_IND}.bootstrap
echo $BS_OUTPUT

echo "generating bootstraps for ${IND}"
${MSMCTOOLS}/multihetsep_bootstrap.py -s 3000000 --out_dir_prefix $BS_OUTPUT --files $BS_INPUT

cd ${MSMCDIR}/bootstrap
ls -d *${POP_OR_IND}.bootstrap_* > ${MSMCDIR}/bs_file_lists/${POP_OR_IND}.bs_file_list.txt

echo "Finished generating bootstraps for ${POP_OR_IND}"

# Run MSMC on all Bootstrapped datasets for your given pop/individual by iteratively sending batch jobs with msmc_4_bootstraps.sh
echo "Running MSMC on newly generated Bootstraps for ${POP_OR_IND}"
for boot in `cat ${MSMCDIR}/bs_file_lists/${POP_OR_IND}.bs_file_list.txt`; do
    echo "running msmc2 on bootstraps for $boot"
    
	sbatch --account=mcnew \
	--job-name=msmc_run.${boot} \
    --partition=standard \
	--output=boot_outs/stdout_${boot} \
    --error=boot_outs/stderr_${boot} \
	--nodes=1 \
	--ntasks=${THREADS} \
	--time=48:00:00 \
	${SCRIPTDIR}/Genomics-Main/B_Phylogenetics/msmc/msmc_4_run_bootstraps.sh -p ${PARAMS} -m ${MSMCPARAMS} -b ${boot} -i ${POP_OR_IND}
done
