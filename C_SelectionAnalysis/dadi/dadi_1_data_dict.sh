#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_data_dictionary
#SBATCH --account=mcnew
##SBATCH --partition=gpu_standard
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
##SBATCH --gres=gpu:1
#SBATCH -o dd_generation.out
#SBATCH -e dd_generation.err
#SBATCH --constraint=hi_mem
#SBATCH --mem-per-cpu=41gb

# --------------------
### Code Section
# --------------------
echo 'Moving to primary directory...'
cd /xdisk/mcnew/finches/ljvossler/finches/dadi/

source /xdisk/mcnew/finches/ljvossler/finches/dadi/param_files/params_dadi.sh
source /xdisk/mcnew/finches/ljvossler/finches/dadi/param_files/params_base.sh

echo 'Activating dadi Micromamba Environment'
source ~/.bashrc
micromamba activate dadi_env

export OUTDIR
export OUT_FOLDER
export VCF_PATH
export POP_PATH

python3 -c 'import os, pickle, dadi; pickle.dump(dadi.Misc.make_data_dict_vcf(os.environ["VCF_PATH"], os.environ["POP_PATH"], calc_coverage=True), open(os.environ["OUTDIR"] + os.environ["OUT_FOLDER"] + "/dd.pkl", "wb"), 2)'
