#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=generate_bootstraps
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
##SBATCH --gres=gpu:1
#SBATCH -o start_bootstrapping.out
#SBATCH -e start_bootstrapping.err
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

for i in `cat samplecodes.txt`; 
    do echo $i;
          sbatch --account=mcnew \
	  --job-name=bootstrap_$i \
	  --partition=standard \
          --nodes=1 \
          --time=24:00:00 \
          --output=outs/stdout_bootstrap_$i \
          --error=outs/stderr_bootstrap_$i \
          --constraint=hi_mem \
          --mem-per-cpu=32gb \
          msmc2_scripts/msmc_4_bootstraps.sh $i; 
done
