echo "Running GONE for $OUTPREFIX with $NUMIND individuals with recombination rate of $RECOMB_RATE..."
cd GONE2/
./gone2 $INDIR/$OUTPREFIX.ped -g $GENO_DTYPE -r $RECOMB_RATE -i $NUMIND -t 4 -o $OUTPREFIX $EXTRA_FLAGS

echo "Organizing Output Files..."

mv ${OUTPREFIX}_GONE2_d2 ../$OUTDIR/
mv ${OUTPREFIX}_GONE2_Ne ../$OUTDIR/
mv ${OUTPREFIX}_GONE2_STATS ../$OUTDIR/


echo "Completed GONE Analysis for $OUTPREFIX"
 
