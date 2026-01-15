#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=submit_gone2_array
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=1
##SBATCH --gres=gpu:1
#SBATCH --output submit_gone2_array.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb
#SBATCH --array=1-6

ARRAY_NAME=gone_run3 # Be sure to rename this parameter if you wish to keep mulitple gone2 runs. Otherwise, all files in an existing array folder will be overwritten.

INPUTPARAMS="$( sed "${SLURM_ARRAY_TASK_ID}q;d" INPUTPOPS )"
OUTFILE="$( echo gone2_${INPUTPARAMS}.out)"

source $INPUTPARAMS

# Create output directory if it doesn't exist
OUTDIR=$ARRAY_NAME/$OUTPREFIX/gone_output
if [ ! -d "$OUTDIR" ]; then
  echo "Directory for gone_output does not exist. Creating it now..."
  mkdir -p "$OUTDIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone_output already exists."
fi

source run_gone2.sh > $OUTDIR/$OUTFILE