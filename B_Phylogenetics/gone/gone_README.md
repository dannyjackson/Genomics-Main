# GONE Pipeline - Estimating Recent Effective Population Size
===========================================
## Script Requirements

Estimating Recent Effective Population Size (Ne) over time in GONE should follow the outline below. Note that GONE is different from other programs in this repo such that it will want your input files and param-file to be in the same directory as the executables. Review the GONE docs if confused.

===========================================
# Step-By-Step Pipeline

**Generating Input Files**
The necessary GONE input files are .ped and .map PLINK files containing pedigree info for genotype calls and variant info respectively. You can generate these from phased or unphased VCF files using `placeholder_script_name.sh`

`placeholder_script_name.sh` will use vcftools to convert your VCF into PLINK .map and .ped files for a given population. See the documentation in the script for an overview of the needed flags/files/parameters to properly convert VCFs to .ped/.map files that GONE will accept. (If it doesn't work for you, go to GONE documentation to better understand what formats it requires and start tinkering.)


**Running GONE**

Once you've generated your input files, ensure that they are in the same directory as the executable GONE scripts that you pulled from the main GONE Git-Repo.

*NOTE: Recall from the GONE docs that all the needed executables must be given permission to run. You can do this on the HPC using `chmod +x EXECUTABLE_FILE_NAME`*

**Plotting Outputs**

WIP


===========================================
# Additional Notes / Key Resources

- *Repo*: https://github.com/esrud/GONE

- *GONE Paper*: "Recent demographic history inferred by high-resolution analysis of linkage disequilibrium" by Enrique Santiago, Irene Novo, Antonio F. Pardiñas, María Saura, Jinliang Wang and Armando Caballero. Molecular Biology and Evolution, 2020 Volume 37, Issue 12, Pages 3642–3653, https://doi.org/10.1093/molbev/msaa169

- *PLINK File Formats Overview*: https://www.cog-genomics.org/plink/1.9/formats

- *VCFtools Manual*: https://vcftools.sourceforge.net/man_latest.html

===========================================
