#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_gone
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=05:00:00
##SBATCH --gres=gpu:1
#SBATCH -o run_gone.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

source GONE/Linux/script_GONE.sh cra_post
