# Repository Overview 
This repository contains generalized scripts and codes used across my multiple genomics projects.
 - For each project, read every script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd.
 - The params.sh file is essential for all scripts in this repository. It is designed to be editable between projects, such that all downstream analysis will run if this file is up-to-date.
   - **The params.sh file does not change all parameters!** Only major ones that change between projects, such as naming conventions, paths to the reference genome, etc! Read all scripts and edit the options/flags of each command according to the needs of your project!
 - This repository should be used as a template for a submodule within any particular project. To do this, we are going to branch this repository and make local edits, but never push our local edits to the parent directory. Use the following steps:

Add the repository as a submodule
```
git submodule add https://github.com/dannyjackson/Genomics-Main/ path/to/project_repository
git submodule update --init --recursive
```
Use a separate branch for local edits
   - Inside the submodule (project_repository/Genomics-Main), create a new branch for each project where you make your local edits.
```
cd path/to/submodule
git checkout -b cardinalis-edit
```
Keep the submodule up to date with the main repository
   - Periodically pull changes from the main Genomics-Main repository into the submodule (without pushing your local changes).
```
cd path/to/submodule
git checkout main
git pull origin main
```
## Repository Conventions
Each module within the repository contains a markdown file in the structured format: **[Module].md** These markdown files describe how the scripts can be called and provides example slurm scripts used to call them on the University of Arizona HPC. 


Every script name also follows a structured format:
**[Module]\_[Step Number]_[Description]**

   - **Module** describes the statistic/program central to the analysis.
   - The **number** indicates the order in which the scripts should be excecuted within that module.
   - The **description** briefly summarizes the purpose of the script.


## Repository Modules

**F<sub>ST</sub> scripts**  
   - These scripts compute Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. They requires SAF files as input, which can be generated using the <scriptname> scripts. They will compute average genome-wide Fst, sliding window Fst, and fst for each SNP, including plots. 

**dxy scripts**  
   - These scripts compute dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD. They requires SAF files as input, which can be generated using the <scriptname> scripts. They will compute average genome-wide dxy, sliding window dxy, and dxy for each SNP, including plots.

**MSMC scripts (in development)**  
   - These scripts analyze a set of genomes using MSMC. They will generate individual plots, population level plots, estimate divergence between populations, and will assess confidence in these outputs using bootstraps


