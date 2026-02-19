#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_ind_mask_cra_pre
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
##SBATCH --gres=gpu:1
#SBATCH -o ind_mask_cra_pre.out

Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ../params_base_df.sh -i /xdisk/mcnew/finches/ljvossler/finches/referencelists/