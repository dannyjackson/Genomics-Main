'''
Author: <Logan Vossler>

==========================================================
Make 2-Dimensional Demography Model
==========================================================

Description:
This script uses dadi to construct 2D models to predict demographic events between species.
Model parameters, optimizations, and uncertainty analyses will be saved in relevant model results directories

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - Bootstrapped SFS files in proper result directory for use in Godambe Uncertainty Analysis
    - If using LowPass workflow, your data dictionary generated in dadi_make_sfs.py is needed in the dadi_results directory
'''


# Required Modules
#==========================================================
import dadi, nlopt, glob, os, json5, sys
import matplotlib.pyplot as plt
import pickle as pkl
from dadi.LowPass import LowPass
from pathlib import Path
import numpy as np


# Function Definitions
#==========================================================
def iso_inbreeding(params, ns, pts):
    '''
    This function can be used to simulate the diverge of two diploid populations and 
    the presence of inbreeding within them. The function is outlined in DADI documentation
    but is not included in the out-of-the-box models, so we must explicity define it here.
    Parameters:
        params: A 5 element list of arbitrary parameters including (in order):
            T: Time of split
            nu1 & nu 2: Respective sizes of both populations
            F1 & F2: Respective inbreeding coefficients for both populations
        ns: Sample sizes
        pts: Number of grid points for the model when plotted
    Returns:
        fs: A frequency spectrum object
    '''
    T, nu1, nu2, F1, F2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2)
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs

def make_2d_demo_model(fs, pop_ids, dadi_model, model_dir, start_params, l_bounds, u_bounds, num_opt=20, lowpass=False):
    '''
    This function makes a 2D demographic model. After optimizing the model, 
    it will save the model parameters into an output file in the results directory.
    Parameters:
        fs: An SFS object containing data across 2 populations
        pop_ids: A 2 element list containing strings of species names
        dadi_model: A string specifying the desired dadi model to run. This will also be used as model results folder name
        model_dir: A string representing the model_specific dadi results directory
        start_params: List of floats/ints of arbitrary starting parameters for model to optimize
        l_bounds: List of floats/ints of lower param optimization bounds
        u_bounds: List of floats/ints of upper param optimization bounds
        num_opt: Integer specifying the number of model optimizations to perform (Defaults to 20 for testing)
        lowpass: Boolean value specifying if Lowpass workflow is desired (Defaults to False to not use lowpass worflow)
    Returns:
        popt: List of int/floats of Optimal model parameters
        model_fs: SFS object created based on optimized model parameters
        model_ex: Our final model function
        pts: List of integers of number of grid points for model
    '''
    # Get Sample Sizes
    n = fs.sample_sizes
    n_total = np.sum(n)
    # Get number of grid points for the model
    pts = [max(n)+20, max(n)+30, max(n)+40]

    # State our model
    if dadi_model == 'iso_inbreeding':
        model = iso_inbreeding
    else:
        model = getattr(dadi.Demographics2D, dadi_model)

    if not fs.folded:
        # Since we have unfolded data, we will wrap the model in a function that adds a parameter to estimate misidentification rate
        model = dadi.Numerics.make_anc_state_misid_func(model)

    # Add another model wrapper that allows for the use of grid points
    model_ex = dadi.Numerics.make_extrap_func(model)

    if lowpass:
        # Calculate Coverage Distribution from data dictionary
        with open('dadi_results/dd.bpkl', 'rb') as file:
            dd = pkl.load(file)
        cov_dist =  LowPass.compute_cov_dist(dd, fs.pop_ids)
        model_ex = LowPass.make_low_pass_func_GATK_multisample(model_ex, cov_dist, fs.pop_ids, nseq=n_total, nsub=n_total, sim_threshold=1e-2, Fx=None)

    # Create a file to store the fitted parameters in your current working directory
    try:
        output = open(model_dir + '_'.join(pop_ids) +'_fits.txt','a')
    except:
        output = open(model_dir + '_'.join(pop_ids) +'_fits.txt','w')
    
    # This is where we run the optimization for our arbitrary model parameters
    # By the end, we will (presumably) get some optimal parameters and the log-liklihood for how optimal they are
    for i in range(num_opt):
        # Slightly alter parameters
        p0 = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=u_bounds,lower_bound=l_bounds)
        popt, ll_model = dadi.Inference.opt(p0, fs, model_ex, pts, lower_bound=l_bounds, upper_bound=u_bounds,algorithm=nlopt.LN_BOBYQA,maxeval=400, verbose=100)
        # Calculate the synonymous theta
        # Also finding optimal scaling factor for data
        model_fs = model_ex(popt, n, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

    # Write results to output file
    res = [ll_model] + list(popt) + [theta0]
    output.write('\t'.join([str(ele) for ele in res])+'\n')
    output.close()

    return popt, model_fs, model_ex, pts
    
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
    # Perform Godambe Uncertainty Analysis
    boots_fids = glob.glob('dadi_results/bootstraps/' + '_'.join(pop_ids) + '/' + '_'.join(pop_ids) + 'boots*.fs')
    boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

    # Godambe uncertainties
    # Will contain uncertainties for the
    # estimated demographic parameters and theta.

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(pop_ids) +'confidence_intervals.txt','w')
    fi.write('Optimized parameters: {0}\n\n'.format(popt))

    # we want to try a few different step sizes (eps) to see if uncertainties very wildly with changes to step size.
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

def compare_sfs_plots(data_fs, model_fs, pop_ids, model_dir):
    '''
    This function plots a comparison spectra between the data and model.
    Will be useful in visually determining model accuracy.
    Parameters:
        data_fs: The actual allele frequency spectra from our samples
        model_fs: The proposed allele frequency spectra from our samples
        pop_ids: A 2 element list containing strings for species being plotted. yaxis is the first list element. xaxis is the second.
        folder_name: string representing the name of folder to put model outputs
    Returns:
        None
    '''
    comp_plot = dadi.Plotting.plot_2d_comp_multinom(model_fs, data_fs, pop_ids=pop_ids)
    plt.savefig(model_dir + '_'.join(pop_ids) + '_comp_plot.png')
    plt.clf()


# Main
#==========================================================
def main():
    '''
    1) Make data dictionary from VCF file
    2) Make the spectrum objects for each species comparison from data dictionary
    3) Plot SFS and Save plots to files
    4) Save SFS objects to files
    5) Make bootstraps
    6) Make demography model SFS and save fits to files
    7) Make SFS / Model SFS comparison plots
    8) Save model SFS to files
    9) Perform uncertainty analysis on all models and save them to files
    '''
    #========================================
    # Import base parameters from user-inputted params_base.sh file
    base_params = Path(sys.argv[1]).read_text().strip().split('\n')
    for line in base_params:
        if 'OUTDIR=' in line: outdir = line.split('=')[1].split('#')[0].strip()

    # Import dadi-specific parameters for model creation from user-inputted dadi_params.json file
    with open(sys.argv[2], 'r') as file:
        dadi_params = json5.load(file)
    dadi_model = dadi_params['DADI MODEL']
    model_params = dadi_params['MODEL PARAMS']
    godambe_eps = dadi_params['GODAMBE STEP SIZES']
    num_opt = dadi_params['PARAM OPTIMIZATIONS']
    lowpass = dadi_params['LOWPASS']

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    result_dir = outdir + 'dadi_results/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    model_dir = result_dir + dadi_model + '/'
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    #========================================

    for dct in model_params:
        pop_ids = model_params[dct][0]
        data_fs = model_params[dct][1]
        start_params = model_params[dct][2]
        l_bounds = model_params[dct][3]
        u_bounds = model_params[dct][4]

        # Make model SFS objects for each species comparison
        popt, model_fs, model_ex, pts = make_2d_demo_model(data_fs, pop_ids, dadi_model, model_dir, start_params, l_bounds, u_bounds, num_opt, lowpass)

        # Plot SFS model/data comparison
        compare_sfs_plots(data_fs, model_fs, pop_ids, model_dir)

        # Save model SFS to files
        model_fs.to_file(model_dir + '_'.join(pop_ids) + 'model_fs')

        # Perform dadi Godambe Uncertainty Analysis
        godambe(popt, pop_ids, model_ex, pts, data_fs, model_dir, godambe_eps)


if __name__ == '__main__':
    main()
