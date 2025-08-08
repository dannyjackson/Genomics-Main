#!/bin/bash
# --------------------
### Directives Section
# --------------------
#SBATCH --job-name=filter_vcf_by_depth
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00

# --------------------
### Code Section
# --------------------
module load bcftools/1.19
module load vcftools/0.1.16

vcftools --vcf "/xdisk/mcnew/finches/ljvossler/finches/dadi_genom_main/vcfs/post_qualitysort.vcf" --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out "/xdisk/mcnew/finches/ljvossler/finches/dadi_genom_main/vcfs/post_depth_filter"
