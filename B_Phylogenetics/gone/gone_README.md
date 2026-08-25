# GONE2 Pipeline - Estimating Recent Effective Population Size
===========================================

Estimating Recent Effective Population Size (Ne) over time in GONE2 should follow the outline below.

## Script Requirements
- Run `gone_setup.sh` prior to performing any input generation or analyses. NOTE: Some GONE2 scripts may require executable permissions to run. `chmod +x file_name`
- Depending on sample size, you may need to boost the memory limit/max loci allowed by GONE2. See the developer docs and `gone_setup.sh` for examples.

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

Run `gone_1_input_generation.sh` to generate a chromosome-subsetted VCF file and plink .ped counterpart. Chromosome subsetting is almost always required since GONE2 only works on scaffolds at least 20cM in length. This is a basic script that does no additional filtering for missing data, MAF, or snp subsetting (all of which may be required for your data. Customize as needed.)

The .map file generated in this step will not contain info on genetic distances (cM) per SNP. If you do not have a suitable reference linkage map to manually interpolate cM values for your SNPs, then this map is useless and you'll need to use a constant recombination rate in GONE2 analysis.

**Running GONE2**

Once you've generated your input files, you can run GONE using `gone_2_analysis.sh`. This script can accept either your .vcf or .ped file inputs depending on type of genotype data you specify.

**Plotting Outputs**

Ne estimates are stored under the filename `outname_GONE2_Ne` file, a tab-delimited file denoting the mean estimates of population size for each generation. You can use `gone_3_plot.r` for a basic example of how to plot this output.

Additional run stats are stored in `popname_GONE2_STATS`.


===========================================
# Additional Notes

- Original GONE Publication: https://doi.org/10.1093/molbev/msaa169
- GONE1 Repo: https://github.com/esrud/GONE
- GONE2 Repo: https://github.com/esrud/GONE2.git

===========================================
