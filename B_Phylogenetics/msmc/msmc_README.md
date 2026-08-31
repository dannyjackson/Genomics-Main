# MSMC2 Pipeline - Estimating Effective Population Size
===========================================

Estimating Effective Population Size (Ne) over time in MSMC2 should follow the outline below.

## Script Requirements
- Run `msmc_setup.sh` to set up directories and create msmc environment. 
- `params_msmc.sh` is present and sources the main `params_base.sh` file.
- These scripts use the pre-compiled version of MSMC2 for Linux, downloadable using micromamba

===========================================
# Step-By-Step Pipeline

**Input File Generation**
You will need:

1) A reference genome mappability mask (`A2.3_generate_mask.sh`).

2) A set of mask and VCF files per chromosome for each individual (if using whatshap phasing pipeline, you can generate these files using `A2.4_individual_mask_vcf.sh` and phase using whatshap `A2.7_phasing_whatshap.sh`.

4) File containing list of sample codes for your population. 

5) Generate MSMC inputs with multi or single individuals using `msmc_2_generateinput_multiInd.sh` or `msmc_2_generateInput_singleInd.sh` respectively. (Phasing required for multi-individual runs)


**MSMC Analysis**

Run MSMC using `msmc_3_runMSMC.sh`. This script supports both single and multi-individual runs. The CPU count for more efficient runs can roughly match your organism's chromosome count.

*NOTE on Haplotype Indices:* For single individual runs (on diploid organisms), only use two indices (usually 0,1). For multi-individual runs, MSMC2 can technically run on as many haplotypes as you want; however, the computational time increases drastically with more haplotypes. 8 haplotypes (4 individuals) is fairly standard. Between 12-16 haplotypes (6-8 diploid individuals) is widely considered the maximum reasonable amount with runtime considerations. Therefore, if you have more samples than this, you will need to select a subset of haplotypes to run.


**MSMC Bootstrapping**

`msmc_4_generate_bootstraps.sh` will create 20 bootstrapped sets of input files for a given individual or population. It will then automatically submit a slurm job and call on `msmc_4_run_bootstraps.sh` to start running MSMC on each bootstrapped input in a batch array. Note that ``msmc_4_run_bootstraps.sh` can be easily run on its own (useful if you already have generated bootstrapped sets and don't wish to waste resources regenerating them).


**Plotting Outputs**

Reference the example script, `msmc_5_plotmsmc.r` for ways to plot output.


===========================================
# Additional Notes

- MSMC Publication: https://doi.org/10.1007/978-1-0716-0199-0_20
- MSMC Repo: https://github.com/stschiff/msmc2
- MSMC-TOOLS Repo: https://github.com/stschiff/msmc-tools/tree/master

- This workflow is a modified form of Jessi Rick's pipeline (https://github.com/jessicarick/msmc2_scripts/).

===========================================
