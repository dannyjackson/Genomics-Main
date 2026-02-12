# GONE2 Pipeline - Estimating Recent Effective Population Size
===========================================

Estimating Recent Effective Population Size (Ne) over time in GONE2 should follow the outline below.

## Script Requirements
- Run `gone2_setup.sh` prior to performing any input generation or analyses. NOTE: Some GONE2 scripts may require executable permissions to run. To do this you may need to run `chmod r+x file_name`
- Properly formatted .ped and .map input files (GONE2 also permits VCFs as inputs if this better fits your workflow)

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

Although GONE2 can accept a number of input filetypes, this pipeline is built for PLINK .ped/.map files. This can be generated from VCFs using `input_generation.sh`. 

 Depending on your scaffold lengths, not all scaffolds may be long enough to use, so you'll likely need to trim a list down to scaffolds greater than 20 centimorgans (about 20 million bp)

NOTE: This input generation script may not be "one-size-fits-all" depending on your data. Be sure to review the `PLINK` and `VCFtools` docs and adjust script parameters as needed to get GONE2 to work with your data.

**Running GONE2**

Once you've generated your input files, you can run GONE using `run_gone2.sh`. When running on the UA HPC, the Puma cluster is usually required since other clusters do not have necessary libraries.
Be sure to revise your script parameters prior to running. Refer to GONE2 repo tutorial docs to understand how to best adjust params based on your needs.


**Plotting Outputs**

Ne estimates are stored under the filename `popname_GONE2_Ne` file. This file contains tab separated columns denoting the mean estimates of population size for each generation. If following the full pipeline, they will be located in your `analyses/gone2_outputs` directory. You can use `plot_gone2.r` for a basic example of how to plot this output.

Some useful information about the analysis is stored in `popname_GONE2_STATS`.


===========================================
# Additional Notes

- Original GONE Publication: https://doi.org/10.1093/molbev/msaa169
- GONE1 Repo: https://github.com/esrud/GONE
- GONE2 Repo: https://github.com/esrud/GONE2.git
- For help with building SLURM arrays on UA HPC: https://hpcdocs.hpc.arizona.edu/running_jobs/batch_jobs/array_jobs/

- If you do not have much genotype and recombination data in cM present in your .map file you'll probably need to set a constant recombination rate in GONE2.
- It is recommended to run these scripts in a batch array, especially when performing analyses across many populations. Use `template_pop_params` for example on how to pass and update parameters for each population in the array.

===========================================
