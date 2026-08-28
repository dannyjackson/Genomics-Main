# MSMC2 Pipeline - Estimating Effective Population Size
===========================================

Estimating Effective Population Size (Ne) over time in MSMC2 should follow the outline below.

## Script Requirements
- Run `msmc_setup.sh` prior to performing any input generation or analyses to set up directories and create your msmc environment. 
- `params_msmc.sh` is present and sources the main `params_base.sh` file.
- These scripts use the pre-compiled version of MSMC2 for Linux, downloadable using micromamba

===========================================
# Step-By-Step Pipeline

**Generating Input Files**
Generating the MSMC haplotype input files requires some preprocessing of your data. You will need:

1) A reference genome mappability mask (generated using `A2.3_generate_mask.sh`).

2) A set of mask and VCF files per chromosome for each individual (generated using `A2.4_individual_mask_vcf.sh`). These VCF files must be phased prior to generating MSMC inputs. Can phase using whatshap (`A2.5_phasing.sh`) or another method of your choice.

4) Create a file (`POP_IND.txt`) containing a list of sample codes for each population you are analyzing (as seen in `sample_POP_IND.txt`). 

Generate MSMC inputs with multi or single individuals using `msmc_2_generateinput_multiInd.sh` or `msmc_2_generateInput_singleInd.sh` respectively. (It is recommended when testing your pipeline to run on single individual first.)


**Running MSMC**

Run MSMC using `msmc_3_runMSMC.sh`. This script supports both single and multi-individual runs. The CPU count for more efficient runs can roughly match your organism's chromosome count.

*NOTE on Haplotype Indices:* For single individual runs (on diploid organisms), only use two indices (usually 0,1). For multi-individual runs, MSMC is designed for up to 12 haplotypes (6 diploid individuals) and cannot handle more than this. Therefore, if you have more than 6 individuals in your population, you will need to select a subset of 12 haplotypes to run.


**Generating Bootstrap Outputs**

`msmc_4_generate_bootstraps.sh` will create 20 bootstrapped sets of input files for a given individual or population. It will then automatically submit a slurm job and call on `msmc_4_run_bootstraps.sh` to start running MSMC on each bootstrapped input in a batch array. Note that ``msmc_4_run_bootstraps.sh` can be easily run on its own (useful if you already have generated bootstrapped sets and don't wish to waste resources regenerating them).


**Plotting Outputs**

An example script, `msmc_5_plotmsmc.r` can be used to plot your outputs locally in RStudio.


===========================================
# Additional Notes

- MSMC Publication: https://doi.org/10.1007/978-1-0716-0199-0_20
- MSMC Repo: https://github.com/stschiff/msmc2
- MSMC-TOOLS Repo: https://github.com/stschiff/msmc-tools/tree/master

- This workflow is a modified form of Jessi Rick's pipeline (https://github.com/jessicarick/msmc2_scripts/).

===========================================
