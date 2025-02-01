# Repository Overview 
This repository contains generalized scripts and codes used across multiple genomics projects.
 - For each project, carefully read every script and revise as needed!
 - The base_setup.sh and params_*.sh files is essential for all scripts in this repository. Each project is associated with a single base_params.sh file that contains parameters that apply across all analyses (e.g. /path/to/reference/genome) and a base_setup.sh file that creates directories and reference lists that are relevant across all modules.
 - Each module additionally has a unique params_*.sh file that contains parameters, options, and flags required for that specific module (i.e. params_msmc.sh).
 - This repository should be used as a template for a submodule within any particular project. To do this, we are going to branch this repository and make local edits, but never push our local edits to the parent directory. Use the following steps:

Add the repository as a submodule:
```
git submodule add git@github.com:dannyjackson/Genomics-Main.git
git submodule update --init --recursive
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

**dxy scripts (in development phase)**  
   - These scripts compute dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD. They requires SAF files as input, which can be generated using the <scriptname> scripts. They will compute average genome-wide dxy, sliding window dxy, and dxy for each SNP, including plots.

**MSMC scripts (in testing phase)**  
   - These scripts analyze a set of genomes using MSMC. They will generate individual plots, population level plots, estimate divergence between populations, and will assess confidence in these outputs using bootstraps


### Contributing to this respository

You do not have to (and should not) edit these files to complete a genomic analysis! Local edits that differ from the Genomics-Main repository will not be reflected in the repository release that we create when we submit for publication. Changes when implementing these modules should only be made to parameter and base_setup files that have been copied to the project's repository.

If you find an error OR if you develop scripts for a new analysis that can be added to this repository, these edits should be made directly to the Genomics-Main repository. This way, they will be subsequently pulled into all projects that rely on this repository. To make and submit edits, do the following:

Ensure that your local submodule is up to date with the latest version of the main Genomics-Main repository:
```
cd path/to/submodule
git checkout main
git pull origin main
```
Inside the submodule (project_repository/Genomics-Main), create a new branch for each project where you make your local edits.
```
cd path/to/submodule
git checkout -b <descriptive name of branch> # e.g. feature-new-analysis or fix-bug-filename
```
Make changes using your preferred text editor and save the files. Then commit any changes that you've made to your branch.
```
git branch # confirms that you're in the branch that you've created
git add .  # Stages all changes
git commit -m "Added new analysis function"
```
Push the branch to GitHub.
```
git push --set-upstream origin feature-new-analysis # if this is the first push you've made to the branch
git push origin feature-new-analysis # if you've pushed the branch before
```