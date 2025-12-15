# GONE2 Pipeline - Estimating Recent Effective Population Size
===========================================

Estimating Recent Effective Population Size (Ne) over time in GONE2 should follow the outline below.

## Script Requirements
- The main workhorse scripts that run GONE are located in the `GONE2` submodule from https://github.com/esrud/GONE2
- Ensure that all scripts located under `./GONE2/` subdirectory are permitted to be executable. To do this you may need to run `chmod r+x file_name`
- Properly formatted .ped and .map input files (GONE2 also permits VCFs as inputs if this better fits your workflow)

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

GONE requires plink `.ped` and `.map` variant files that can be generated using `input_file_generation.sh`. This script requires a VCF file for your population and a file containing the scaffolds you wish to include in analyses separated by line (as seen in `example_scaffold_sets.txt`). Depending on your scaffold lengths, not all scaffolds may be long enough to use, so you'll likely need to trim this list down to scaffolds greater than 20 centimorgans (about 20 million bp)

NOTE: This input generation script may not be "one-size-fits-all" depending on your data. Be sure to review the `plink` and `vcftools`docs and adjust script parameters as needed.

**Running GONE2**

Once you've generated your input files, you can run GONE using `submit_gone2.sh` to send a slurm job to the HPC's Puma cluster. (If another cluster works, then that's fine, but Puma seems to usually have all the libraries required to run GONE2)
Be sure to revise your script parameters prior to running. Refer to GONE2 repo tutorial docs to understand how to best adjust params based on your needs.


**Plotting Outputs**

GONE outputs can be found in `popname_GONE2_Ne` file. This file contains tab separated columns denoting the geometric mean estimates of population size for each generation.

You can use `plot_gone.r` to plot and compare time-series outputs


===========================================
# Additional Notes

- Original GONE Publication: https://doi.org/10.1093/molbev/msaa169
- GONE1 Repo: https://github.com/esrud/GONE

===========================================
