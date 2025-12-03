#!/bin/sh

# This script uses VCFtools and PLINK to generate properly formatted .ped and .map input files for use in GONE analysis of recent effective population size
# Since this process is usually fairly quick, this script is designed to just be run in an interactive session, but edit as needed

module load vcftools
module load plink

VCFFILE=$1
SCAFFOLDS=$2
OUTPREFIX=$3

if [ ! -d "$OUTPREFIX/gone_input" ]; then
  echo "Directory for gone_input does not exist. Creating it now..."
  mkdir -p "$OUTPREFIX/gone_input" # -p creates parent directories if they don't exist
else
  echo "Directory for gone_input already exists."
fi

# Convert VCF to proper .map file format for GONE
vcftools --gzvcf $VCFFILE --plink --chrom-map $SCAFFOLDS --remove-indels --out $OUTPREFIX
rm $OUTPREFIX.ped
mv $OUTPREFIX.map $OUTPREFIX/gone_input

# Convert VCF to proper .binary file format
plink --vcf $VCFFILE --make-bed --allow-extra-chr --snps-only --out $OUTPREFIX
# Convert binary file format to text .ped file format
plink --bfile $OUTPREFIX --recode --allow-extra-chr --out $OUTPREFIX
mv $OUTPREFIX.ped $OUTPREFIX/gone_input

rm $OUTPREFIX.map

mv $OUTPREFIX.bed $OUTPREFIX
mv $OUTPREFIX.log $OUTPREFIX
mv $OUTPREFIX.nosex $OUTPREFIX
mv $OUTPREFIX.bim $OUTPREFIX
mv $OUTPREFIX.fam $OUTPREFIX


