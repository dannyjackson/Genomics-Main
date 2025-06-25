#!/bin/bash

#echo 'Activating dadi micromamba environment...'
source ~/.bashrc
micromamba activate dadi_env

echo '=========================='
echo JOB_NAME $1
echo NUM_OPT $2
echo LOWPASS $3
echo PLOT DEMES $4
echo OUTDIR $5
echo CURRENT MODEL $6
echo POP IDS "${7} ${8}"
echo SFS PATH $9
echo $MODEL_JSON ${10}
echo '=========================='

source /xdisk/mcnew/finches/ljvossler/finches/dadi/param_files/params_dadi.sh

START_PARAMS=$(python -c "import json; print(json.load(open('${10}'))['${6}'][0])")
LOW_PARAMS=$(python -c "import json; print(json.load(open('${10}'))['${6}'][1])")
UP_PARAMS=$(python -c "import json; print(json.load(open('${10}'))['${6}'][2])")

echo $GIM_STEPS

echo $START_PARAMS
echo $LOW_PARAMS
echo $UP_PARAMS
echo '=========================='

echo 'Starting dadi script...'
python3 /xdisk/mcnew/finches/ljvossler/finches/dadi/scripts/Genomics-Main/C_SelectionAnalysis/dadi/model_making_test.py -j $1 -f ${11} -n $2 -l -d -o $5 -m $6 -p "${7} ${8}" -s $9 --start_params "${START_PARAMS}" --low_params "${LOW_PARAMS}" --up_params "${UP_PARAMS}" --eps "$GIM_STEPS"