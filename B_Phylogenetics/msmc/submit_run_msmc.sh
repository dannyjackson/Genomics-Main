#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_msmc_pop1
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=24:00:00
##SBATCH --gres=gpu:1
#SBATCH -o run_msmc_pop1.out
#SBATCH -e run_msmc_pop1.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

PARAMS=/path/to/params_base.sh

source ${PARAMS}

source ${SCRIPTDIR}/B_Phylogenetics/msmc/msmc_3_runMSMC.sh -p "${PARAMS}" -m params_msmc.sh -i pop1
