# Important dadi Script Information

## Script Requirements
Each dadi script requires the base_params.sh file and dadi_params.json file.
If you need additional explanation on the basic JSON file structure, read the json_file_structure_overview.md.

Each script will require the following python modules (and any additional dependencies for them)
*All other unstated modules should be installed by default on most Python installations*

- dadi (Check Official Documentation to ensure all dependencies met)
- json5 (This is an alternative Python JSON file parser that allows for comments in JSON file structure) <-- If you run into issues with this module, edit the scripts to use the default json parser instead (import json) and remove all comments from dadi_params.json

When executing in an interactive terminal or submitting as a bash script to an HPC, each params file must be submitted as arguments when executing a dadi script.
*Use Python Version 3.9 for smoother experience.*
```
module load python/3.9/3.9.10

python3 dadi_make_sfs.py base_params.sh dadi_params.json
```

## Pipeline
Making Demography models in dadi should follow the outline below.
 1) Generate and filter your VCF as needed. We have found that dadi prefers bcftools for this step over other programs such as ANGSD
 2) Run dadi_make_sfs.py to create your Site Frequency Spectra.
 3) Run dadi_make_2dmodels.py to generate your demography models.
 4) Run dadi_run_godambe.py to perform dadi GIM Uncertainty Analysis


## Parameter Input Notes
Below are notes and recommendations on how some parameters should be stated to properly run dadi scripts:

**File Path Parameters**:
All file path parameters should be strings representing either the full or relative file path.

**2D Site Frequency Spectra**:
2D-Site Frequency Spectra parameters should contain your population IDs ("strings"), your population projections (integers), and a T/F boolean statement of whether the SFS should be polarized. All SFS params should be placed under the key, SFS PARAMS, in the event of creating multiple SFS's.
Example of an SFS Param List:
```
"SFS_1":[["pop0_ID", "pop1_ID"], [22, 33], false]
```

**Bootstrapped SFS Params**:
These parameters should be integers stating the number of bootstrapped datasets you want to generate, and the chunk size of each one (in number of base pairs).
Example:
```
"BOOTSTRAP PARAMS": [100, 1e7]
```

**dadi Model Selection**:
Your desired dadi model MUST be a 2D model offered in dadi.Demographics2D submodule. State model as actual method name (Refer to dadi documentation)
Example: If your desired model function is dadi.Demographics2D.split_mig, then state just the *split_mig* portion.
```
"DADI MODEL": "split_mig"
```
**Note, dadi_make_2dmodel.py has a custom an iso_inbreeding() model function that is not offered in dadi-out-of-the-box function. If you wish to use it, state "iso_inbreeding" as the parameter**

**LowPass Workflow**:
If your data is 10x coverage or lower, you may plan to use dadi's LowPass workflow. If so, change the boolean parameter to true.

**Model Parameter Optimizations**:
Recommended to run 20 rounds of optimizations for testing and 100 rounds on HPC for full analyses.
Number of runs should be stated as an integer

**2D Demography Models**:
2D-Models should contain population IDs, SFS file paths, and lists of your starting, lower, and bounds of model parameters. Your pop_ids and SFS file paths should be strings, your parameters should be lists of floats or integers.
All Model params should be placed under the key, MODEL PARAMS, in the event of running multiple models.
Recall that if your SFS is unfolded, you will need one additional value at the end of your param lists that acts as the misidentification parameter.
Example: 
```
"Model_1":[["pop0_ID", "pop1_ID"], "path/to/fs", [22, 33, 44, 55], [0.01, 0.01, 0.01, 0.01], [100, 100, 100, 100]]
```
    
**Godambe Uncertainty Analysis**:
You can change the range of your step sizes that you test depending on your data. (Refer to dadi documentation and discussion posts)
Steps sizes should be stated as lists of floats.



