#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=start_multiple_demo_models
#SBATCH --account=mcnew
##SBATCH --partition=gpu_standard
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
##SBATCH --gres=gpu:1
#SBATCH -o start_model_runs.out
#SBATCH -e start_model_runs.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

# --------------------
### Code Section
# --------------------
echo 'Moving to primary directory...'
cd /xdisk/mcnew/finches/ljvossler/finches/dadi/

echo 'Load Param Files'
source param_files/params_base.sh
source param_files/params_dadi.sh

echo 'Loading Python 3.11.4...'
module load python/3.11/3.11.4

echo $POP_IDS

for i in `seq 0 $(python -c "import json; print(len(json.load(open('${MODEL_JSON}')).keys())-1)")`; 
    do echo $i
    CURRENT_MODEL=$(python -c "import json; print(list(json.load(open('${MODEL_JSON}')).keys())[$i])")
    echo '============================'
    echo 'Current Model: ' $CURRENT_MODEL
    sbatch --account=mcnew \
    --job-name=${JOB_NAME}.${CURRENT_MODEL} \
    --partition=standard \
    --output=scripts/dadi_model_outs/dadi_model.${JOB_NAME}.${CURRENT_MODEL}.out \
    --error=scripts/dadi_model_outs/dadi_model.${JOB_NAME}.${CURRENT_MODEL}.err \
    --nodes=1 \
    --ntasks=1 \
    --time=04:00:00 \
    /xdisk/mcnew/finches/ljvossler/finches/dadi/scripts/2d_demo_model.sh ${JOB_NAME} ${NUM_OPT} ${LOWPASS} ${PLOT_DEMES} ${OUTDIR} ${CURRENT_MODEL} ${POP_IDS} ${SFS_PATH} ${MODEL_JSON} ${OUT_FOLDER}
    echo '============================'
done