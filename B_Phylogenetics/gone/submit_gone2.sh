#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_gone2
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=4
##SBATCH --gres=gpu:1
#SBATCH --output run_gone2.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

#=PARAMS===========================
POPNAME=pop1
NUMIND=5
POP_PATH=../$POPNAME
RECOMB_RATE=3.1
OUTDIR=$POP_PATH/gone_output
#==================================

echo "Running GONE for $POPNAME with $NUMIND individuals with recombination rate of $RECOMB_RATE..."
cd GONE2/
./gone2 $POP_PATH/gone_input/$POPNAME.ped -g 0 -r $RECOMB_RATE -i $NUMIND -t 4 -o $POPNAME


echo "Organizing Output Files..."

if [ ! -d $OUTDIR ]; then
  echo "Output directory for $POPNAME does not exist. Creating it now..."
  mkdir -p $OUTDIR
else
  echo "Output directory for $POPNAME already exists. Moving on..."
fi
mv ${POPNAME}_GONE2_d2 $OUTDIR/
mv ${POPNAME}_GONE2_Ne $OUTDIR/
mv ${POPNAME}_GONE2_STATS $OUTDIR/


echo "Completed GONE Analysis for $POPNAME"
 
