# dadi Script Requirements
Each dadi script requires the base_params.sh file and dadi_params.json file.

When executing in an interactive terminal or submitting as a bash script to an HPC, each params file must be submitted as arguments when executing a dadi script.

```
module load python/3.11.4

python3 dadi_make_sfs.py base_params.sh dadi_params.json
```

## dadi Parameter File Structure
The dadi-specific parameters for making SFS and/or 2D demography models are contained in a JSON file. The basic structure is as follows:
```
{
    KEY_1: Value_1,
    KEY_2: Value_2,
    KEY_3: Value_3
}
```
**NOTE:** dadi_params.json has comments in it denoted by double-slashes ('//').
The dadi scripts use the alternative Python JSON parser (json5) that is able to correctly read these comments. However, if you encounter errors with reading comments, then remove all comments before running the scripts.