#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=run_gone_cra_post
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
##SBATCH --gres=gpu:1
#SBATCH --output run_gone_cra_post.out
##SBATCH --constraint=hi_mem
##SBATCH --mem-per-cpu=41gb

POPNAME=cra_post

GONEDIR=GONE/Linux


echo "Copying Datafiles..."
# Copying needed datafiles to required GONE directory
cp $POPNAME/gone_input/$POPNAME.ped $GONEDIR
cp $POPNAME/gone_input/$POPNAME.map $GONEDIR
cp params_gone.sh $GONEDIR/INPUT_PARAMETERS_FILE

echo "Running GONE for $POPNAME..."
cd $GONEDIR
source script_GONE.sh $POPNAME


echo "Organizing Output Files..."

rm $POPNAME.map
rm $POPNAME.ped
rm -r TEMPORARY_FILES

cd ../../
mkdir $POPNAME/gone_output
mv $GONEDIR/outfileHWD $POPNAME/gone_output/
mv $GONEDIR/Output_d2_$POPNAME $POPNAME/gone_output/
mv $GONEDIR/Output_Ne_$POPNAME $POPNAME/gone_output/
mv $GONEDIR/OUTPUT_$POPNAME $POPNAME/gone_output/
mv $GONEDIR/timefile $POPNAME/gone_output/
mv $GONEDIR/seedfile $POPNAME/gone_output/

echo "Completed GONE Analysis for $POPNAME"
 
