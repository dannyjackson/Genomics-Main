#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_input_multi_cra_pre
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
##SBATCH --gres=gpu:1
#SBATCH -o input_multi_cra_pre.out

Genomics-Main/B_Phylogenetics/msmc/msmc_2_generateInput_multiInd.sh -p params_msmc.sh -i /xdisk/mcnew/finches/ljvossler/finches/referencelists/CRA_pre_IND.txt -s CRA_pre