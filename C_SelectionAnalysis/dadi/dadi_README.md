# dadi Script Requirements
Each dadi script requires the base_params.sh file and dadi_params.json file.
If you need additional explanation on the basic JSON file structure, read the json_file_structure_overview.md

When executing in an interactive terminal or submitting as a bash script to an HPC, each params file must be submitted as arguments when executing a dadi script.

```
module load python/3.11.4

python3 dadi_make_sfs.py base_params.sh dadi_params.json
```

**NOTE:** dadi_params.json has comments in it denoted by double-slashes ('//'). Comments are not supported in basic Python JSON Parsers. To get around this, the dadi scripts use the alternative Python JSON parser (json5) that is able to correctly read these comments. However, if you still encounter errors with reading comments, them before running the scripts.

## Parameter Input Notes

**2D Site Frequency Spectra**
2D-Site Frequency Spectra parameters should contain your population IDs, your population projections, and a T/F boolean statement of whether the SFS should be polarized.
Example:
```
"SFS_1":[["pop0_ID", "pop1_ID"], [22, 33], false]
```

**dadi Model Selection**
Desired dadi model MUST be a 2D model offered in dadi.Demographics2D submodule. State model as actual method name (Refer to dadi documentation)
Example: If your desired model function is dadi.Demographics2D.split_mig, then state just the *split_mig* portion.
```
"DADI MODEL": "split_mig"
```

**Model Parameter Optimizations**
Recommended to run 20 rounds of optimizations for testing and 100 rounds on HPC for full analyses

**2D Demography Models**
2D-Models should contain population IDs, SFS paths, and lists of your starting, lower, and bounds of model parameters
Recall that if your SFS is unfolded, you will need one additional value at the end of your param lists that acts as the misidentification parameter.
    Example: 
```
"Model_1":[["pop0_ID", "pop1_ID"], "path/to/fs", [22, 33, 44, 55], [0.01, 0.01, 0.01, 0.01], [100, 100, 100, 100]]
```
    
**Godambe Uncertainty Analysis**
You can change the range of your step sizes that you test depending on your data. (Refer to dadi documentation and discussion posts)



