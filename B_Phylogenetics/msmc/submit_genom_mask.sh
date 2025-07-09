#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=make_genom_mask
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=26:00:00
##SBATCH --gres=gpu:1
#SBATCH -o make_mask.out
#SBATCH -e make_mask.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb


source ~/.bashrc
micromamba activate mask_env


# Make file location edits and such directly within the .py script in this case. It takes no parameters
python msmc-tools/makeMappabilityMask.py
