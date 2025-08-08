#!/bin/bash

echo 'Activating dadi micromamba environment...'
source ~/.bashrc
micromamba activate dadi_env

echo '=========================='
echo JOB NAME $1
echo NUM OPT $2
echo LOWPASS $3
echo PLOT DEMES $4
echo OUTDIR $5
echo CURRENT MODEL $6
echo SFS PATH $7
echo MODEL JSON ${8}
echo OUTFOLDER ${9}
echo ALIGNMENT LEN ${10}
echo TOTAL SNPS ${11}
echo MU ${12}
echo GENERATION TIME ${13}
echo '=========================='

source param_files/params_dadi.sh

START_PARAMS=$(python -c "import json; print(json.load(open('${8}'))['${6}'][0])")
LOW_PARAMS=$(python -c "import json; print(json.load(open('${8}'))['${6}'][1])")
UP_PARAMS=$(python -c "import json; print(json.load(open('${8}'))['${6}'][2])")

echo $GIM_STEPS

echo $START_PARAMS
echo $LOW_PARAMS
echo $UP_PARAMS
echo '=========================='

echo 'Starting dadi script...'
python3 dadi_3_demo_model.py -j $1 -f ${9} -n $2 -l -d -o $5 -m $6 -s $7 --start_params "${START_PARAMS}" --low_params "${LOW_PARAMS}" --up_params "${UP_PARAMS}" --eps "$GIM_STEPS" --align_len ${10} --num_snps ${11} --mu ${12} -g ${13}