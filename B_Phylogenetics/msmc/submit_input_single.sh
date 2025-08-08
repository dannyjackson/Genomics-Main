#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_input_single
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=26:00:00
##SBATCH --gres=gpu:1
#SBATCH -o input_single.out
#SBATCH -e input_single.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

source ~/.bashrc
micromamba activate msmc_env

source params_msmc.sh

source msmc2_scripts/msmc_3_generateinput_singleindv.sh -p params_base.sh
