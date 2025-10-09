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

PARAMS=/path/to/params_base.sh

source ${PARAMS}

source ${SCRIPTDIR}/B_Phylogenetics/msmc/msmc_2_generateInput_singleInd.sh -p "${PARAMS}" -m params_msmc.sh
