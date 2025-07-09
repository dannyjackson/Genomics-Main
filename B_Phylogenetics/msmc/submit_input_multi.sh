#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_input_multi
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
##SBATCH --gres=gpu:1
#SBATCH -o input_multi.out
#SBATCH -e input_multi.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

for POP in `cat popcodes.txt`;
	do echo $POP
	sbatch --account=mcnew \
	--job-name=msmc_run.${POP} \
    --partition=standard \
	--output=slurm_output/msmc_input.${POP}.%j \
	--nodes=1 \
	--ntasks=1 \
	--time=10:00:00 \
    msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo ${POP}_IND.txt` `echo $POP`
done

