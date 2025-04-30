'''
Author: <Logan Vossler>

==========================================================
Make 2-Dimensional Demography Model
==========================================================

Description:
This script uses dadi to construct 2D models to predict demographic events between species.
We may want to edit it for 1D functionality, but not yet. 
Model parameters, optimizations, and uncertainty analyses will be saved in relevant model results directories

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - Bootstrapped SFS files in proper result directory for use in Godambe Uncertainty Analysis
    - If using LowPass workflow, your data dictionary generated in dadi_make_sfs.py is needed in the dadi_results directory
    - If parallelizing this script using CUDA Enabled GPUs (recommended for if doing many runs), you will also need Python's PYCUDA 
        and scikit-cuda modules (along with loading the Nvidia CUDA toolkits from an HPC)
'''


# Required Modules
#==========================================================
import dadi, nlopt, os, sys
import matplotlib.pyplot as plt
import dill as pkl
from pathlib import Path
from dadi.LowPass import LowPass # Can comment out if not using LowPass workflow
import json5 # Can switch to normal json module if this one causes issues
import demes
import demesdraw # Can comment out if not making deme plots


# Function Definitions
#==========================================================

# Model Functions
#========================
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


# Model Optimization
#========================
def make_2d_demo_model(fs, pop_ids, dadi_model, model_dir, result_dir, start_params, l_bounds, u_bounds, num_opt=20, lowpass=False):
    '''
    This function makes a 2D demographic model. After optimizing the model, 
    it will save the model parameters into an output file in the results directory.
    Parameters:
        fs: An SFS object containing data across 2 populations
        pop_ids: A 2 element list containing strings of species names
        dadi_model: A string specifying the desired dadi model to run. This will also be used as model results folder name
        model_dir: A string representing the model_specific dadi results directory
        result_dir: A string representing the dadi results directory
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
    print('---> Make Model Function')
    model_ex = dadi.Numerics.make_extrap_func(model)

    if lowpass:
        print('---> Diverting to LowPass workflow...')
        print('---> Loading Coverage Distribution...')
        with open(result_dir + '_'.join(pop_ids) + '_cov_dist.pkl', 'rb') as file:
            cov_dist = pkl.load(file)
        print('---> Making LowPass Model Function...')
        model_ex = LowPass.make_low_pass_func_GATK_multisample(model_ex, cov_dist, pop_ids, nseq=n, nsub=n, sim_threshold=1e-2, Fx=None)

    # Create a file to store the fitted parameters in your current working directory
    try:
        output = open(model_dir + '_'.join(pop_ids) + '_fits.txt','a')
    except:
        output = open(model_dir + '_'.join(pop_ids) + '_fits.txt','w')
    
    # This is where we run the optimization for our arbitrary model parameters
    # By the end, we will (presumably) get some optimal parameters and the log-liklihood for how optimal they are
    print('---> Optimizing Starting Parameters...')
    for i in range(num_opt):
        # Slightly alter parameters
        p0 = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=u_bounds,lower_bound=l_bounds)
        popt, ll_model = dadi.Inference.opt(p0, fs, model_ex, pts, lower_bound=l_bounds, upper_bound=u_bounds, algorithm=nlopt.LN_BOBYQA, maxeval=400, verbose=100)
        # Calculate the synonymous theta
        # Also finding optimal scaling factor for data
        model_fs = model_ex(popt, n, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

    # Write results to output file
    print('---> Saving Optimized Model Parameters...')
    res = [ll_model] + list(popt) + [theta0]
    output.write('\t'.join([str(ele) for ele in res])+'\n')
    output.close()

    return popt, model_fs, model_ex, pts


# Plotting Functions
#========================
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

def make_demes_plot(model_dir, pop_ids):
    '''
    This function makes a demes plot from dadi model information
    Parameters:
        model_dir: A string representing the model_specific dadi results directory
        pop_ids: A 2 element list containing strings of species names
    Returns:
        None
    '''
    deme_model = dadi.Demes.output(Nref=200)
    deme_plot = demesdraw.tubes(deme_model)
    deme_plot.figure.savefig(model_dir + '_'.join(pop_ids) + '_demes_plot.png')
    plt.clf()

# Main
#==========================================================
def main():
    '''
    1) Import Base Parameters
    2) Import Dadi-SFS Specific Parameters
    3) Check for required directories
    6) Generate demography model SFS and save fits to files
    7) Plot SFS / Model SFS comparisons/residuals
    8) Save model SFS to file
    9) Save uncertainty analysis on all models to file
    '''
    #========================================
    # Basic check for enough parameter files inputted
    print('Checking File Arguments...')
    if len(sys.argv) != 3:
        raise IndexError('Not enough arguments. Did you include both the Base and dadi-specific parameter files?')
    
    # Import base parameters from user-inputted params_base.sh file
    print('Storing Needed Base Parameters...')
    base_params = Path(sys.argv[1]).read_text().strip().split('\n')
    outdir = None
    for line in base_params:
        if 'OUTDIR=' in line: outdir = line.split('=')[1].split('#')[0].strip() 
    if not outdir:
        raise NameError("Couldn't find Out-Directory. Check that OUTDIR is specified as in example Genomics-Main base_params.sh.")

    # Import dadi-specific parameters for model creation from user-inputted dadi_params.json file
    print('Storing Needed dadi Parameters...')
    with open(sys.argv[2], 'r') as file:
        dadi_params = json5.load(file)
    job_name = dadi_params['JOB NAME']
    dadi_model = dadi_params['DADI MODEL']
    model_params = dadi_params['MODEL PARAMS']
    num_opt = dadi_params['PARAM OPTIMIZATIONS']
    lowpass = dadi_params['LOWPASS']
    # Check if CUDA Enabled
    dadi.cuda_enabled(dadi_params['CUDA ENABLED'])
    print('dadi.cuda status: ' + str(dadi.cuda_enabled()))

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    print('Verifying Directories...')
    # If using lowpass, make a lowpass directory
    result_dir = outdir + job_name + '/lowpass/' if lowpass else outdir + 'dadi_results/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    model_dir = result_dir + dadi_model + '/'
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    
    #========================================
    # Enter a loop to perform model-making runs for each species combo
    gim_params = []
    for dct in model_params:
        print('\nGetting ' + dct + ' Parameters...')
        pop_ids, fs_fname, start_params, l_bounds, u_bounds = model_params[dct]
        data_fs = dadi.Spectrum.from_file(fs_fname)

        # Make model SFS objects for each species comparison
        print('Generating ' + dct + '...')
        popt, model_fs, model_ex, pts = make_2d_demo_model(data_fs, pop_ids, dadi_model, model_dir, result_dir, start_params, l_bounds, u_bounds, num_opt, lowpass)
        gim_params.append([popt, pop_ids, model_ex, pts, data_fs])

        # Plot SFS model/data comparison
        print('Plotting SFS Comparison for ' + dct + '...')
        compare_sfs_plots(data_fs, model_fs, pop_ids, model_dir)

        # Save model SFS to files
        print('Saving SFS file for ' + dct + '...')
        model_fs.to_file(model_dir + '_'.join(pop_ids) + '_model_fs')
    
    # Save GIM Params to .pkl file for Uncertainty Analysis
    print('\nSaving Intermediate GIM Params file to Model Directory...')
    with open(model_dir + 'gim_params.pkl', 'wb') as file:
        pkl.dump(gim_params, file)

    print('\n**Model Creation Complete**')

if __name__ == '__main__':
    main()
