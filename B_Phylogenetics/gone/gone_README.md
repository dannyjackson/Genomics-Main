# GONE2 Pipeline - Estimating Recent Effective Population Size
===========================================

Estimating Recent Effective Population Size (Ne) over time in GONE2 should follow the outline below.

## Script Requirements
- The GONE2 executables should be located in `GONE2` submodule from https://github.com/esrud/GONE2. Ensure that programs here are permitted to be executable. To do this you may need to run `chmod r+x file_name`
- Properly formatted .ped and .map input files (GONE2 also permits VCFs as inputs if this better fits your workflow)

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

GONE2 takes a number of possible input files, most commonly .ped/.map and VCF files. .ped/.map variant files that can be generated from VCFs using `input_generation.sh`. A file containing the scaffolds you wish to include in analyses separated by line (as seen in `example_scaffold_sets.txt`) can be used to grab a subset of chromosomes. Depending on your scaffold lengths, not all scaffolds may be long enough to use, so you'll likely need to trim this list down to scaffolds greater than 20 centimorgans (about 20 million bp)

NOTE: This input generation script may not be "one-size-fits-all" depending on your data. Be sure to review the `plink` and `vcftools`docs and adjust script parameters as needed to get GONE2 to work with your data.

**Running GONE2**

Once you've generated your input files, you can run GONE using `run_gone2.sh`. Note that when running on the UA HPC, the Puma cluster is usually required since other clusters do not have necessary libraries for GONE2.
Be sure to revise your script parameters prior to running. Refer to GONE2 repo tutorial docs to understand how to best adjust params based on your needs.


**Plotting Outputs**

GONE outputs can be found in `popname_GONE2_Ne` file. This file contains tab separated columns denoting the geometric mean estimates of population size for each generation.
You can use `plot_gone2.r` to plot and compare time-series outputs. NOTE: Depending on the flags used during GONE2 analysis, the output filename structure can sometimes change. In this case, `run_gone2.sh` will not detect or move the output-Ne file into the organized directory, but you can find it in the GONE2/ submodule directory.


===========================================
# Additional Notes

- Original GONE Publication: https://doi.org/10.1093/molbev/msaa169
- GONE1 Repo: https://github.com/esrud/GONE

- If you do not have much genotype and recombination data in cM present in your .map file you'll probably need to set a constant recombination rate in GONE2.
- It is recommended to run these scripts in a batch array, especially when performing anaylses across many populations.

===========================================
