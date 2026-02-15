#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_lrt_analysis
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00

# --------------------
### Code Section
# --------------------
echo 'Moving to primary directory...'
cd /xdisk/mcnew/finches/ljvossler/finches/dadi/

echo 'Load Param Files'
source param_files/params_base.sh
source param_files/params_dadi.sh

echo 'Activating dadi Micromamba Environment'
source ~/.bashrc
micromamba activate dadi_env

python3 dadi_3_LRT.py -j ${JOB_NAME} -s ${SFS_PATH} -o ${OUTDIR} -b ${BOOT_DIR}  -i "${NESTED_INDICES}" -n ${NULL_MODEL} -t ${TEST_MODEL} --null_popt "${NULL_POPT}" --test_popt "${TEST_POPT}"
