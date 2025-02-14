## dadi Parameter File Structure
The dadi-specific parameters for making SFS and/or 2D demography models are contained in a JSON file. The basic structure is as follows:

All JSON data is placed within curly brackets {} and is formatted in key:value pairs. The values can be user-edited, but the keys should **NEVER** be changed by the user.

## Basic Formatting
- All text string values must be placed in double quotes. 
- No quotes are needed to state float/integer values. 
- true/false statements are Boolean values and do not require double quotes (ensure that they are all lowercase).
- We can assign lists of values to a single key with square brackets []
```
{
    "VCF PATH": "path/to/vcf-file"
    "PARAM OPTIMIZATIONS": 20,
    "LOWPASS": false,
    "GODAMBE STEP SIZES":[0.1, 0.01, 0.001, 0.0001, 0.00001]
}
```
## Nested JSON Structure
- We have nested JSON structures to store SFS and model parameter data. See below how we have a key (SFS PARAMS) that stores more key value pairs that contain our parameters for creating each SFS that we want.
- Note that we can use lists to assign different value types to a single key. See how each of our SFS keys contain nested lists holding strings, integers, and boolean values for use in creating an SFS.
```
{
    "SFS PARAMS": {
        "SFS_1":[["pop0_ID", "pop1_ID"], ["pop0_projection", "pop1_projection"], false],
        "SFS_2":[["pop0_ID", "pop1_ID"], ["pop0_projection", "pop1_projection"], false],
        "SFS_3":[["pop0_ID", "pop1_ID"], ["pop0_projection", "pop1_projection"], false]
                  }
}
```