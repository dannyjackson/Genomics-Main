'''
Author: <Logan Vossler>

==========================================================
Run Likelihood-Ratio Test
==========================================================

Description:
-----
'''


# Required Modules
#==========================================================
import dadi, glob, os, sys
import dill as pkl
from pathlib import Path
import json5 # Can switch to normal json module if this one causes issues


# Function Definitions
#==========================================================
def load_bootstraps(result_dir, pop_ids):
    '''
    This is a helper function that gets the bootstrap data for our uncertainty analyses.
    Parameters:
        result_dir: A string representing the dadi results directory
        pop_ids: A list containing strings of population names
    Returns:
        boots_syn: A list of bootstrapped SFS objects
    '''
    # Get Bootstrapped datasets
    boots_fids = glob.glob(result_dir + 'bootstraps/' + '_'.join(pop_ids) + '/' + '_'.join(pop_ids) + 'boots*.fs')
    boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]
    return boots_syn

def likelihood(popt, model_ex, pts, fs, model_dir, eps, result_dir, pop_ids, lrt_indices):
    '''
    This function performs dadi likelihood ratio test on 2D demographic models.
    It requires SFS bootstraps generated from dadi_make_sfs.py and will save the results to text files in the dadi model results directory.
    Parameters:
        popt: List of int/floats of Optimal model parameters
        model_ex: Our final model function
        pts: ist of integers of number of grid points for model
        fs: An SFS object containing data across 2 populations
        model_dir: A string representing the model_specific dadi results directory
        eps: List of floats of step sizes to test in our confidence intervals
        result_dir: A string representing the dadi results directory
        pop_ids: A 2 element list containing strings of species names
        lrt_indices: List of indices to fix for simple model
    Returns:
        None
    '''
    # Get Bootstrapped datasets
    boots_syn = load_bootstraps(result_dir, pop_ids)

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(pop_ids) +'LRT_results.txt','w')

    for steps in eps:
        adj = dadi.Godambe.LRT_adjust(func_ex=model_ex, grid_pts=pts, all_boot=boots_syn, p0=popt, data=fs, nested_indices=lrt_indices, multinom = True, eps=steps)
    fi.close()

# Main
#==========================================================
def main():
    '''
    1) Import Base Parameters
    2) Import Dadi-SFS Specific Parameters
    3) Check for required directories
    4) Load GIM Parameters
    5) Perform GIM Uncertainty Analysis and write results to files
    '''
    #========================================
    # Basic check for enough parameter files inputted
    print('Checking File Arguments...')
    if len(sys.argv) != 3:
        raise IndexError('Not enough arguments. Did you include the Base Params and dadi-specific params files?')
    
    # Import base parameters from user-inputted params_base.sh file
    print('Storing Needed Base Parameters...')
    base_params = Path(sys.argv[1]).read_text().strip().split('\n')
    outdir = None
    for line in base_params:
        if 'OUTDIR=' in line: outdir = line.split('=')[1].split('#')[0].strip() 
    if not outdir:
        raise NameError("Couldn't find Out-Directory. Check that OUTDIR is specified as in example Genomics-Main base_params.sh.")

    # Import dadi-specific parameters needed for GIM Analysis from user-inputted dadi_params.json file
    print('Storing Needed dadi Parameters...')
    with open(sys.argv[2], 'r') as file:
        dadi_params = json5.load(file)
    job_name = dadi_params['JOB NAME']
    dadi_model = dadi_params['DADI MODEL']
    lowpass = dadi_params['LOWPASS']

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    print('Verifying Directories...')
    # If using lowpass, make a lowpass directory inside specified results folder
    result_dir = outdir + job_name + '/lowpass/' if lowpass else outdir + job_name + '/'
    if not os.path.exists(result_dir):
        raise FileNotFoundError("Couldn't find specified result directory. Check that specified Result Directory matches the previously generated directory in data SFS creation.")

    # Specify model directory
    model_dir = result_dir + dadi_model + '/'
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)