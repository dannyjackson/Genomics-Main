# WIP LDHelmet Pipeline - Estimating Fine-Scale Recombination Rates
===========================================

Estimating fine-scale recombination rates in LDHelmet should follow the outline below.

## Script Requirements
- Run `ldhelmet_setup.sh` prior to performing any input generation or analyses to set up directories and create your LDHelmet environment. This will require some manual work. Refer to the script for more info.
- `params_ldhelmet.sh` is present and references the main `params_base.sh` file.

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

The primary inputs accepted by LDHelmet are fasta files, the creation and format of which is explained more by the developers. This pipeline instead uses `vcftools` to convert phased VCF files to the alternative `.snp` and `.loc` inputs. This can be done using `ldhelmet_input_generation.sh`

**Running LDHelmet**

You can run LDHelmet using `ldhelmet_analysis.sh`. Be sure to revise your `params_ldhelmet.sh` script parameters prior to running. Refer to the LDHelmet repo docs to understand how to best adjust params based on your needs.


**Result Post-Processing**

LDHelmet will output a binary file with the extension `.post`. You can textualize and assess these outputs using `ldhelmet_postprocess.sh`. Note that this step often requires a moderate-to-high amount of memory.


===========================================
# Additional Notes

- LDHelmet Publication: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003090
- LDHelmet Repo: https://github.com/popgenmethods/LDhelmet
- Original LDHat Repo: https://github.com/auton1/LDhat/tree/master
- For help with building SLURM arrays on UA HPC: https://hpcdocs.hpc.arizona.edu/running_jobs/batch_jobs/array_jobs/


- It is recommended to run these scripts in a batch array, especially when performing analyses across many populations.

===========================================
