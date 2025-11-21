# GONE Pipeline - Estimating Recent Effective Population Size
===========================================
## Script Requirements
- Properly formatted .ped and .map files present in `gone/` directory.
- Script parameters stated in `params_gone.sh`

The main workhorse scripts that run GONE are located in the `GONE/` submodule. This module is a forked version of GONE. It is functionally identical, with only minor edits to allow for running in a SLURM job on the UA HPC.

Forked Repo: 

Original Repo: https://github.com/esrud/GONE

Estimating Recent Effective Population Size (Ne) over time in GONE should follow the outline below. All scripts will require edits to a copy of the params_gone.sh file to function

===========================================
# Step-By-Step Pipeline

**Generating Input Files**

GONE requires plink `.ped` and `.map` variant files in specific formats that can be generated using `input_file_generation.sh`. This script requires a VCF file for your population and a Scaffold Map file that maps your scaffolds to numbers (as seen in `example_scaffolds.txt`)

**Running GONE**

Once you've generated your input files, you can run GONE using `submit_gone.sh` to send a slurm job to the HPC.
Be sure to revise your script parameters in `params_gone.sh` prior to running.


**Plotting Outputs**

GONE Ne outputs can be found in `Output_Ne_popname` file. This file contains tab separated columns denoted the geometric mean estimates of population size for each generation.

You can use `gone_plot.r` to plot your outputs locally in RStudio. (The GONE output files are small and it can useful to tweak plot parameters on the fly in RStudio)


===========================================
# Additional Notes

-

===========================================
