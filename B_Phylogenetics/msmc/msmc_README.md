# MSMC2 Pipeline - Estimating Effective Population Size
===========================================
## Script Requirements
Each script requires the base_params.sh file and params_msmc.sh file.

Estimating Effective Population Size (Ne) over time in MSMC2 should follow the outline below. All scripts will require edits to a copy of the params_msmc.sh file to function

===========================================
# Step-By-Step Pipeline

**Generating Input Files**
Generating the MSMC haplotype input files requires some preprocessing of your data. You will need:

1) A reference genome mappability mask (generated using `A2.3_generate_mask.sh`).

2) A set of mask and VCF files per chromosome for each individual (generated using `A2.4_individual_mask_vcf.sh`). 

3) Lastly, these VCF files must be phased (generated using `A2.5_phasing.sh`).

Finally you can choose to generate your input files using either `msmc_2_generateinput_multiInd.sh` or `msmc_2_generateInput_singleInd.sh` depending on if you wish to estimate `Ne` based on single or multi individual haplotype data. (It is recommended when testing your pipeline to run on single individual first.)


**Running MSMC**
Once you've generated your input files, you can run MSMC using `msmc_3_runMSMC.sh`. This script supports both single and multi-individual runs. When doing so, keep in mind that MSMC can run across multiple CPUs. Be sure to edit the `$THREADS` parameter as needed. Note that the amount of CPUs needed for an optimally efficient run usually matches the amount of chromosomes your organism has.
**NOTE:** Lines 104-106 of `msmc_3_runMSMC.sh` are using to define the haplotype indices being analyzed. It is hard-coded to input all possible haplotype pairs into MSMC. Edit a copy of the the script if you wish to change this.

**Generating Bootstrap Outputs**
`msmc_4_generate_bootstraps.sh` will create 20 bootstrapped sets of input files for a given individual or population. It will then call `msmc_4_run_bootstraps.sh` to start running MSMC on each bootstrapped input in separate batch jobs. Note that ``msmc_4_run_bootstraps.sh` can be easily edited to run on its own (useful if you already have generated bootstrapped sets and don't wish to waste resources regenerating them).

**Plotting Outputs**
Use `msmc_5_plotmsmc.r` to plot your outputs locally in RStudio. (The MSMC output files are small and it can useful to tweak plot parameters on the fly in RStudio; however, a version of the plotting script that generates some predefined plots may be created in the future)


===========================================
# Additional Notes

**NOTE:** This workflow is a modified form of Jessi Rick's pipeline (found here: https://github.com/jessicarick/msmc2_scripts/). Her documentation is very well done and can be an additional resource.

===========================================


# TO-DO:

-Integrate input generation scripts better into GenMain workflow
-R Plotting should allow for running on HPC
