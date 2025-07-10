#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=make_sfs_plot
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --constraint=hi_mem
#SBATCH --mem-per-cpu=42gb
##SBATCH --mem=100gb
#SBATCH -o sfs_test.out
#SBATCH -e sfs_test.err

# --------------------
### Code Section
# --------------------
cd /xdisk/mcnew/finches/ljvossler/finches/dadi/

echo 'Activating dadi Micromamba Environment'
source ~/.bashrc
micromamba activate dadi_env

source param_files/params_base.sh
source param_files/params_dadi.sh

echo 'SFS Generation Parameters:'
echo "Job Name: ${JOB_NAME}"
echo "Out Directory: ${OUTDIR}"
echo "Results Folder: ${OUT_FOLDER}"
echo "Pop IDs: ${POP_IDS}"
echo "Number of Chromosomes: ${NUM_CHROMS}"
echo "VCF FilePath: ${VCF_PATH}"
echo "Pop Assignment FilePath: ${POP_PATH}"
echo "Lowpass: ${LOWPASS}"
echo "Polarized: ${POLARIZE}"
echo "Bootstrap Parameters: ${BOOTSTRAP_PARAMS}"

python3 /xdisk/mcnew/finches/ljvossler/finches/dadi/scripts/Genomics-Main/C_SelectionAnalysis/dadi/dadi_2_sfs.py -f ${OUT_FOLDER} -p "${POP_IDS}" -n "${NUM_CHROMS}" -l ${LOWPASS} -o ${OUTDIR} -v ${VCF_PATH} -i ${POP_PATH} -b "${BOOTSTRAP_PARAMS}"

