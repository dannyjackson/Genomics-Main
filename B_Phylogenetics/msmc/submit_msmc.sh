#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_msmc_indv
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=3:00:00
##SBATCH --gres=gpu:1
#SBATCH -o msmc_all_indv.out
#SBATCH -e msmc_all_indv.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

source params_base.sh

for i in `cat ${FILENAME_LIST}`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=msmc_run.${i} \
        --partition=standard \
	--output=slurm_output/msmc_run.${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=10:00:00 \
	msmc2_scripts/msmc_3_runMSMC.sh $i
done
