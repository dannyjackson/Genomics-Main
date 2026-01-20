#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=submit_input_array
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --ntasks-per-node=1
##SBATCH --gres=gpu:1
#SBATCH --output submit_input_array.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb
#SBATCH --array=1-6

ARRAY_NAME=input_run3 # Be sure to rename this parameter if you wish to keep mulitple input runs. Otherwise, all files in an existing array folder will be overwritten.

INPUTPARAMS="$( sed "${SLURM_ARRAY_TASK_ID}q;d" INPUTPOPS )"
OUTFILE="$( echo inputs_${INPUTPARAMS}.out)"

source $INPUTPARAMS

# Create output directory if it doesn't exist
OUTDIR=$ARRAY_NAME/$OUTPREFIX/gone_input
if [ ! -d "$OUTDIR" ]; then
  echo "Directory for gone_input does not exist. Creating it now..."
  mkdir -p "$OUTDIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone_input already exists."
fi

source input_generation.sh > $OUTDIR/$OUTFILE