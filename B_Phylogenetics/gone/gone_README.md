# GONE Pipeline - Estimating Recent Effective Population Size
===========================================

Estimating Recent Effective Population Size (Ne) over time in GONE should follow the outline below. All scripts will require edits to a copy of the params_gone.sh file to function

## Script Requirements
- The main workhorse scripts that run GONE are located in the `GONE/` submodule from https://github.com/esrud/GONE
- Ensure that all scripts located under `./GONE/Linux/PROGRAMMES/` subdirectory are  are permitted to be executable using `chmod r+x file_name`
- Properly formatted .ped and .map files present in `gone/` directory.
- Script parameters stated in `params_gone.sh`

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

GONE requires plink `.ped` and `.map` variant files in specific formats that can be generated using `input_file_generation.sh`. This script requires a VCF file for your population and a Scaffold Map file that maps your scaffolds to integers in tab-separated format (as seen in `example_scaffold_map.txt`)

NOTE: Especially if working with non-model organisms with chromosome names that these software don't immediately recognize, you may need to include flags such as `--allow-extra-chr` in `plink` to avoid errors in reading chromosome names.

**Running GONE**

Once you've generated your input files, you can run GONE using `submit_gone.sh` to send a slurm job to the HPC.
Be sure to revise your script parameters in `params_gone.sh` prior to running.


**Plotting Outputs**

GONE Ne outputs can be found in `Output_Ne_popname` file. This file contains tab separated columns denoting the geometric mean estimates of population size for each generation.

You can use `gone_plot.r` to plot your outputs locally in RStudio.


===========================================
# Additional Notes

- Original GONE Publication: https://doi.org/10.1093/molbev/msaa169

===========================================
