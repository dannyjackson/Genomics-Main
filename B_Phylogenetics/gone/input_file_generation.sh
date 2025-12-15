module load vcftools
module load plink
module load python/3.11.4

VCFFILE=$1
SCAFFOLD_LIST=$2
OUTPREFIX=$3

if [ ! -d "$OUTPREFIX/gone_input" ]; then
  echo "Directory for gone_input does not exist. Creating it now..."
  mkdir -p "$OUTPREFIX/gone_input" # -p creates parent directories if they don't exist
else
  echo "Directory for gone_input already exists."
fi

CHR_SUBSET_FLAGS=$(for name in $(cat $SCAFFOLD_LIST); do echo --chr $name; done)

# Create a filtered VCF to only include a specific subset of chromosomes
vcftools $CHR_SUBSET_FLAGS --gzvcf $VCFFILE --recode --recode-INFO-all --out $OUTPREFIX
#vcftools --chr NC_044571.1 --chr NC_044572.1 --chr NC_044573.1 --chr NC_044574.1 --chr NC_044575.1 --chr NC_044576.1 --chr NC_044577.1 --chr NC_044578.1 --chr NC_044579.1 --chr NC_044580.1 --chr NC_044581.1 --chr NC_044582.1 --gzvcf par_pre.phased.sorted.vcf.gz --recode --recode-INFO-all --out par_pre_subset

# Convert VCF to plink formats
# --allow-extra-chr to deal with non-standard chromosome names
# --thin-count to randomly sample a subset of SNPs from the VCF (needed since GONE2 can't handle too many SNPs by default)
plink --vcf $OUTPREFIX.recode.vcf --thin-count 2000000 --allow-extra-chr --recode --out $OUTPREFIX
# Note that if you don't have info on SNP position in a genetic map (cM), you'll probably need to set a fixed recombination rate when running GONE2 or do some extra analyses to find this info yourself. (In this case the .map file is not useful)

# Due to not having standard chromosome names and that outputted files aren't always consistent with what GONE2 wants, we use these python scripts to reformat the data.
python map_clean.py --map $OUTPREFIX.map
python ped_clean.py --ped $OUTPREFIX.ped

# Organize output files
mv $OUTPREFIX.map $OUTPREFIX/gone_input
mv $OUTPREFIX.ped $OUTPREFIX/gone_input
mv $OUTPREFIX.log $OUTPREFIX/gone_input
mv $OUTPREFIX.nosex $OUTPREFIX/gone_input
mv $OUTPREFIX.recode.vcf $OUTPREFIX/gone_input
