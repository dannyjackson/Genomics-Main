#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_gone2_pipetest
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=4
##SBATCH --gres=gpu:1
#SBATCH --output run_gone2_pipetest.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

POPNAME=pipetest
POP_PATH=../$POPNAME

echo "Running GONE for $POPNAME..."
cd GONE2/
./gone2 $POP_PATH/gone_input/$POPNAME.ped -g 2 -r 1 -t 4 -o $POPNAME


echo "Organizing Output Files..."

if [ ! -d $POP_PATH/gone_output ]; then
  echo "Output directory for $POPNAME does not exist. Creating it now..."
  mkdir -p $POP_PATH/gone_output
else
  echo "Output directory for $POPNAME already exists. Moving on..."
fi
mv $POP_PATH/gone_input/$POPNAME_GONE2_d2 $POP_PATH/gone_output/
mv $POP_PATH/gone_input/$POPNAME_GONE2_Ne $POP_PATH/gone_output/
mv $POP_PATH/gone_input/$POPNAME_GONE2_STATS $POP_PATH/gone_output/


echo "Completed GONE Analysis for $POPNAME"
 
