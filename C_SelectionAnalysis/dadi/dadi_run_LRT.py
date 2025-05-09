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
import dfinbredmodels
import json5 # Can switch to normal json module if this one causes issues


# Function Definitions
#==========================================================
def load_bootstraps(boot_dir, pop_ids):
    '''
    This is a helper function that gets the bootstrap data for our uncertainty analyses.
    Parameters:
        result_dir: A string representing the dadi results directory
        pop_ids: A list containing strings of population names
    Returns:
        boots_syn: A list of bootstrapped SFS objects
    '''
    # Get Bootstrapped datasets
    boots_fids = glob.glob(os.path.join(boot_dir, '*'))
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
    print('Storing Needed dadi_LRT Parameters...')
    with open(sys.argv[2], 'r') as file:
        dadi_LRT_params = json5.load(file)
    data_fs = dadi_LRT_params['DATA SFS']
    test_model = dadi_LRT_params['TEST MODEL'][0].split('.')[-1]
    test_model_params = dadi_LRT_params['TEST MODEL'][1]
    null_model = dadi_LRT_params['NULL MODEL'][0].split('.')[-1]
    null_model_params = dadi_LRT_params['NULL MODEL'][1]
    nested_idx = dadi_LRT_params['NESTED INDICES']
    result_dir = dadi_LRT_params['RESULT OUTDIR']
    boot_dir = dadi_LRT_params['BOOTSTRAP DIR']

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    print('Verifying Directories...')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    #========================================
    # Get function for test model
    try:
        test_func = getattr(dadi.Demographics2D, test_model)
    except:
        test_func = getattr(dfinbredmodels, test_model)
    # Get function for null model
    try:
        null_func = getattr(dadi.Demographics2D, null_model)
    except:
        null_func = getattr(dfinbredmodels, null_model)
    
    # Put wrapper around funcs
    test_ex = dadi.Numerics.make_extrap_log_func(test_func)
    null_ex = dadi.Numerics.make_extrap_log_func(null_func)

    # Get log-likelihoods for each model
    n = data_fs.sample_sizes
    pts = [max(n)+20, max(n)+30, max(n)+40]
    model_test = test_ex(test_model_params, n, pts)
    model_null = null_ex(null_model_params, n, pts)
    ll_test = dadi.Inference.ll_multinom(model_test, data_fs)
    ll_null = dadi.Inference.ll_multinom(model_null, data_fs)

# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple (no-mig) best-fit params.  Since evalution is done with more
# complex model, need to insert zero migration value at corresponding migration
# parameter index in complex model. And we need to tell the LRT adjust function
# that the 3rd parameter (counting from 0) is the nested one.

    # Grab Bootstraps
    boot_syn = load_bootstraps(boot_dir)

    adj = dadi.Godambe.LRT_adjust(test_ex, pts, boot_syn, null_model_params, data_fs, nested_indices=nested_idx, multinom=True)
    D_adj = adj*2*(ll_test - ll_null)
    print('Adjusted D statistic: {0:.4f}'.format(D_adj))

# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
    pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
    print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))