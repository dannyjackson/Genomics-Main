#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
#scriptdir=$(dirname "$0")
#source ${scriptdir}/msmc_params.sh

source ~/.bashrc
micromamba activate msmc_env

source params_msmc.sh
source params_base.sh

IND=$1

# Verify that output location for msmc_outputs exists
if [ -d '${OUTDIR}/bootstrap/outputs/${IND}' ]; then
    echo ${IND} 'bootstrap out directory exists'
else
    mkdir ${OUTDIR}/bootstrap/outputs/${IND}
    mkdir ${OUTDIR}/bootstrap/outputs/${IND}/log_and_loop
    echo ${IND} 'bootstrap out directory created'
fi

#input for the bootstrapping
BS_INPUT=`for s in $(cat SCAFFOLDS.txt); do find ${OUTDIR}/single_indv_data/input_single/ -maxdepth 1 -name "msmc_input.*${IND}.${s}*.txt"; done`

echo 'BS_INPUT======================='
echo $BS_INPUT
echo '======================='

#output from the bootstrapping 
BS_OUTPUT=${OUTDIR}/bootstrap/${IND}.bootstrap
echo $BS_OUTPUT

#echo "generating bootstraps for ${IND}"
msmc-tools/multihetsep_bootstrap.py --out_dir_prefix $BS_OUTPUT --files $BS_INPUT

cd ${OUTDIR}/bootstrap
ls -d *${IND}.bootstrap_* > ${OUTDIR}/bs_file_lists/${IND}.bs_file_list.txt



##### run msmc ####
#MSMC_BS=$(for x in `cat ${OUTDIR}/bs_file_lists/${IND}.bs_file_list.txt`; do find ${OUTDIR}/bootstrap/$x -maxdepth 2 -name "bootstrap_multihetsep*.txt"; done)

for x in `cat ${OUTDIR}/bs_file_lists/${IND}.bs_file_list.txt`; do
    MSMC_BS=$(find ${OUTDIR}/bootstrap/$x -maxdepth 2 -name "bootstrap_multihetsep*.txt")
    echo $MSMC_BS
    
    MSMC_OUTPUT=${OUTDIR}/bootstrap/outputs/${IND}/msmc_output.$x
    echo $MSMC_OUTPUT
    
    echo "running msmc2 on bootstraps for $x"
    msmc2_Linux -o $MSMC_OUTPUT -I 0,1 $MSMC_BS
    
    mv $MSMC_OUTPUT*loop.txt ${OUTDIR}/bootstrap/outputs/${IND}/log_and_loop/
    mv $MSMC_OUTPUT*log ${OUTDIR}/bootstrap/outputs/${IND}/log_and_loop/
    
    echo "done with msmc bootstraps for $x"
    
done
