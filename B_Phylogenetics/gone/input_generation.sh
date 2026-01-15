module load vcftools
module load plink
module load python/3.11.4
module load bcftools

echo ARRAY NAME: $ARRAY_NAME

echo VCF FILE PATH: $VCFFILE
echo SCAFFOLD SUBSET: $SCAFFOLD_LIST
echo POPULATION: $OUTPREFIX
echo MAX MISSING DATA: $MAX_MISSING

OUTDIR=$ARRAY_NAME/$OUTPREFIX/gone_input

if [ ! -d "$OUTDIR" ]; then
  echo "Directory for gone_input does not exist. Creating it now..."
  mkdir -p "$OUTDIR" # -p creates parent directories if they don't exist
else
  echo "Directory for gone_input already exists."
fi

CHR_SUBSET_FLAGS=$(for name in $(cat $SCAFFOLD_LIST); do echo --chr $name; done)

# Create a filtered VCF to only include a specific subset of chromosomes
# --max-missing to set proportion of missing data you'll permit.
vcftools $CHR_SUBSET_FLAGS --gzvcf $VCFFILE --recode --recode-INFO-all --max-missing $MAX_MISSING --out $OUTPREFIX

# Convert VCF to plink formats
# --allow-extra-chr to deal with non-standard chromosome names
# --thin-count to randomly sample a subset of SNPs from the VCF (needed since GONE2 can't handle too many SNPs by default)
plink --vcf $OUTPREFIX.recode.vcf --allow-extra-chr --recode --out $OUTPREFIX
# Note that if you don't have info on SNP position in a genetic map (cM), you'll probably need to set a fixed recombination rate when running GONE2 or do some extra analyses to find this info yourself. (In this case the .map file is not useful)

# Due to not having standard chromosome names and that outputted files aren't always consistent with what GONE2 wants, we use these python scripts to reformat the data. (No filtering or analyses done here, just moving the numbers around)
python map_clean.py --map $OUTPREFIX.map
python ped_clean.py --ped $OUTPREFIX.ped

# Organize output files
mv $OUTPREFIX.map $OUTDIR
mv $OUTPREFIX.ped $OUTDIR
mv $OUTPREFIX.log $OUTDIR
mv $OUTPREFIX.nosex $OUTDIR
mv $OUTPREFIX.recode.vcf $OUTDIR