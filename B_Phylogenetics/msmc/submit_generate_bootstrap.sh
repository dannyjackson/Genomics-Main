#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_bootstrap_pop1
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
##SBATCH --gres=gpu:1
#SBATCH -o start_bootstrapping_pop1.out
#SBATCH -e start_bootstrapping_pop1.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

PARAMS=/path/to/params_base.sh

source ${PARAMS}

source ${SCRIPTDIR}/B_Phylogenetics/msmc/msmc_4_generate_bootstraps.sh -p "${PARAMS}" -m params_msmc.sh -i pop1
