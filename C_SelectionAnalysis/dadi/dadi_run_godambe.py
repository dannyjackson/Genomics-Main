'''
Author: <Logan Vossler>

==========================================================
Run Godambe Uncertainty Analysis
==========================================================

Description:
This script uses dadi to construct 2D models to predict demographic events between species.
Model parameters, optimizations, and uncertainty analyses will be saved in relevant model results directories

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - Bootstrapped SFS files in proper result directory for use in Godambe Uncertainty Analysis
'''


# Required Modules
#==========================================================
import dadi, glob, os, json5, sys
import dill as pkl
from pathlib import Path


# Function Definitions
#==========================================================
def godambe(popt, pop_ids, model_ex, pts, fs, model_dir, eps):
    '''
    This function performs dadi Godambe Uncertainty Analysis on out 2D demographic models.
    It requires SFS bootstraps generated from dadi_make_sfs.py and will save the confidence intervals to text files
    in the dadi model results directory.
    Parameters:
        popt: List of int/floats of Optimal model parameters
        pop_ids: A 2 element list containing strings of species names
        model_ex: Our final model function
        pts: ist of integers of number of grid points for model
        fs: An SFS object containing data across 2 populations
        model_dir: A string representing the model_specific dadi results directory
        eps: List of floats of step sizes to test in our confidence intervals
    Returns:
        None
    '''
    # Get Bootstrapped datasets
    boots_fids = glob.glob('dadi_results/bootstraps/' + '_'.join(pop_ids) + '/' + '_'.join(pop_ids) + 'boots*.fs')
    boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

    # Godambe uncertainties will contain uncertainties for the estimated demographic parameters and theta.

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(pop_ids) +'confidence_intervals.txt','w')
    fi.write('Optimized parameters: {0}\n\n'.format(popt))

    # We want to try a few different step sizes (eps) to see if uncertainties very wildly with changes to step size. (Ideally they shoud not)
    for steps in eps:
        # Get optimzed parameters * 100 (possibly can solve low parameter values leading to floating point arithmetic errors)
        #popt_100 = [param * 100 for param in popt]
        # Do normal uncertainty analysis
        uncerts_adj = dadi.Godambe.GIM_uncert(func_ex=model_ex, grid_pts=pts, all_boot=boots_syn, p0=popt, data=fs, eps=steps, log=True)
        uncerts_str = ',  '.join([str(ele) for ele in uncerts_adj])
        fi.write('Godambe Uncertainty Array Output: [' + uncerts_str + ']\n')
        fi.write('Estimated 95% uncerts (with step size '+str(eps)+'): {0}\n'.format(1.96*uncerts_adj[:-1]))
        fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(popt-1.96*uncerts_adj[:-1]))
        fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(popt+1.96*uncerts_adj[:-1]))
        
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
    dadi_model = dadi_params['DADI MODEL']
    eps = dadi_params['GODAMBE STEP SIZES']

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    print('Verifying Directories...')
    result_dir = outdir + 'dadi_results/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    model_dir = result_dir + dadi_model + '/'
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    #========================================
    # Load GIM Params from intermediate file generated in dadi_make_2d_model.py
    print('\nLoading GIM Parameters from gim_params.pkl...')
    with open(model_dir + 'gim_params.pkl') as file:
        gim_params = pkl.loads(file)
    
    #========================================
    # Enter a loop to perform GIM Analysis for each species combo
    for lst in gim_params:
        popt, pop_ids, model_ex, pts, fs = lst
        print('\nPerforming GIM Analysis for ' + '_'.join(pop_ids) + ' ' + dadi_model +' Model...')
        godambe(popt, pop_ids, model_ex, pts, fs, model_dir, eps)
    print('\n**GIM Analysis Complete**')


if __name__ == '__main__':
    main()