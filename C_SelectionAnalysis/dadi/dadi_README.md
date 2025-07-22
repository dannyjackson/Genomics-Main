# dadi Pipeline - Inferring Demographic Histories
===========================================
## Script Requirements
Each dadi script requires the base_params.sh file and params_dadi.sh file.

## General Pipeline
Making Demography models in dadi should follow the outline below. All scripts will require edits to the params_dadi.sh file to function
 1) Generate and filter your VCF as needed. Any relevant tool should be sufficient for this, but we have found the most success with bcftools over other programs such as ANGSD
 2) Generate and Save VCF data into dadi data dictionaries as .pkl files using dadi_1_data_dict.sh
 3) Create your Site Frequency Spectra using dadi_2_sfs.py.
 3) Optimize your demography models using dadi_3_demo_model.py.
 4) Perform dadi Likelihood Ratio Tests between generated models using dadi_4_LRT.py

===========================================
# Step-By-Step Pipeline

*NOTE: These first two steps can be difficult to run depending on the size of your VCF data and HPC resources available to you. Though the scripts are functional and we have done our best to split the load of data dictionary creation and loading across them, they can commonly run into Out-Of-Memory Errors. Therefore, be prepared to troubleshoot in interactive sessions (and possibly just generate your final data SFS there too, since these sessions have been more consistent than sbatch). On the UofA HPC, both of these steps always require a high memory CPU node*

## Data Dictionary Creation
`dadi_1_data_dict.sh` is used to generate and save VCF data into the python dictionary format the dadi requires. 

## SFS Creation
`dadi_2_sfs.py` is how you will generate your data Site Frequency Spectra and accompanying bootstraps. You can do this for both 1D and 2D SFS. Reminder that you must have previously generated your data dictionary for this script to function. You can run this job on an HPC using `submit_sfs.sh` and editing relevent parameters in `params_dadi.sh`.

## Demography Model Optimization
`dadi_3_demo_model.py` is the workhorse behind creating and optimizing any demography models you run. It can support both 1D and 2D models. It can also generate Demes Plots and will perform dadi's Godambe Uncertainty Analysis (GIM). This GIM analysis is useful for determining if your optimized model parameters are reasonable. However, comparing different models using Likelihood Ratio Tests is highly recommended. 
If you are generating models in bulk, `submit_multi_demo_model.sh` can start multiple model optimization jobs on an HPC. When submitting model jobs in bulk, note the expected use of the model param JSON files (`dadi_1dmodel_params.json` and `dadi_2dmodel_params.json`). Such JSON files are where you can 1) state which models you want to run, and 2) edit their parameter boundaries and starting positions as needed. *NOTE: The model names in JSON must correspond exactly to the function names within either the dadi submodules or any custom modules you create. (For Example: "bottlegrowth_1d" corresponds to the dadi model function `dadi.Demographics1D.bottlegrowth_1d`)*

You can also edit a copy of the `submit_demo_model.sh` script (which is what the former uses to start each job) if you wish to test/rerun a single model at a time.

## Likelihood Ratio Tests
`dadi_4_LRT.py` performs Likelihood Ratio Comparisons between previously generated nested models (comparing a more complex model to a less complex model). You can run this job on an HPC using `submit_LRT.sh` and editing relevent parameters in `params_dadi.sh`. 

NOTE: You cannot just compare any model to another. Rather, they must be "nested", meaning that they are very similar and the complex model represents a special case of the simple model.
For example: You could run an LRT between a bottlegrowth model (simple) and another bottlegrowth model that factors in inbreeding (complex). But you cannot run an LRT between, say, a bottlegrowth model and a split-migration model. Be sure to review the dadi Documenation and Google Group Forum to review how to perform LRTs properly.

===========================================
# Additional Notes
For more info on each of these key points, review the dadi documentation, API, and Google Groups

Docs: https://dadi.readthedocs.io/en/latest/

API: https://dadi.readthedocs.io/en/latest/api/dadi/

Google Groups: https://groups.google.com/g/dadi-user

**LowPass Workflow**:
dadi has a built in submodule (LowPass) to account for biases in model optimization due to low coverage data. If your data is 10x coverage or lower, you may plan to use dadi's LowPass workflow and compare its results to the standard workflow.

**Model Parameter Optimizations**:
dadi developers recommend running at least 20 rounds of optimizations for testing and 100 rounds on HPC for full analyses. For generating final results, going over 100 optimizations won't hurt (other than increasing computational time).

**Demography Models**:
You will likely start out using dadi's builtin 1D or 2D demographic models to get a feel for dadi and your dadi. (The file, dfinbredmodels.py is our example of this for our Darwin Finches project)
If you end up writing some custom model functions, then keep in mind you'll need to edit dadi_3_demo_model.py to import the .py module containing that code.

**GPU Parallelization**:
If you are running on an HPC that has CUDA Enabled Nvidia GPUs, you may wish to speed up model-making in dadi substantially by enabling dadi's CUDA submodule.
Currently these dadi scripts do not implement that feature since we have had trouble getting it working on our HPC and we find that standard CPU model runs with 100 optimizations run in a satisfactory amount of time (ie: A standard 1D dadi model with 100 optimizations can finish within 1 hour).

===========================================

