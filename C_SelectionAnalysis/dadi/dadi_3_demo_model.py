'''
Author: <Logan Vossler>

==========================================================
Optimize Demography Model
==========================================================

Description:
This script uses dadi to construct either 1D or 2D models to predict demographic events between/within species/populations.
Model parameters, optimizations, and uncertainty analyses will be saved in relevant model results directories

Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - Bootstrapped SFS files in proper result directory for use in Godambe Uncertainty Analysis
    - If using LowPass workflow, a coverage distribution generated in dadi_make_sfs.py is needed in the dadi_results directory
    - If generating Demes Plots, the demesdraw module is needed.
    - If parallelizing this script using CUDA Enabled GPUs (recommended for if doing many runs), you will also need Python's PYCUDA 
      and scikit-cuda modules (along with loading the Nvidia CUDA toolkits from an HPC)
      # Note that this functionality has been...finnicky on the HPC. Usually for running the standard dadi models, a normal CPU node is satisfactory in both resources and speed.
'''


# Required Modules
#==========================================================
import argparse
import dadi, nlopt, os, glob
import matplotlib.pyplot as plt
import dfinbredmodels as custom_models # Our .py module containing custom demo-models
import json
#import pycuda
#import skcuda


# CMD-LINE Arguments
#==========================================================
parser=argparse.ArgumentParser()
parser.add_argument("-j", "--job_name", type=str, help='Name of Job')
parser.add_argument("-f", "--out_folder", type=str, help='Name folder where to dump results')
parser.add_argument("-n", "--num_opt", type=int, help='Number of Model Optimizations', default=20)
parser.add_argument("-l", "--lowpass", help='Do Lowpass Pipeline', action='store_true')
parser.add_argument("-d", "--demes", help='Generate Demes Plot', action='store_true')
parser.add_argument("-s", "--sfs_path", type=str, help='Path to Data SFS')
parser.add_argument("--up_params", type=str, help='List of Upper Param Bounds', nargs='+')
parser.add_argument("--low_params", type=str, help='List of Lower Param Bounds', nargs='+')
parser.add_argument("--start_params", type=str, help='List of Starting Params', nargs='+')
parser.add_argument("-o","--outdir", type=str, help='Out Directory Path')
parser.add_argument("-m","--model", type=str, help='Model Function')
parser.add_argument("-e", "--eps", type=str, help='List of step sizes for GIM Uncertainty Analysis', nargs='+')
args = parser.parse_args()


# Function Definitions
#==========================================================
def plot_model_fs(data_fs, model_fs, model_dir, dim=2):
    '''
    Plots a comparison spectra between the data and model. Useful in visually determining model accuracy.
    Parameters:
        data_fs: Data allele frequency spectrum
        model_fs: Proposed Model Allele frequency spectrum
        model_dir: String representing the path to directory of model results
        dim: Int representing dimensionality of sfs (defaults to 2)
    Returns:
        None
    '''
    if dim == 2:
        comp_plot = dadi.Plotting.plot_2d_comp_multinom(model_fs, data_fs, vmin=1, pop_ids=data_fs.pop_ids)
    elif dim == 1:
        comp_plot = dadi.Plotting.plot_1d_comp_multinom(model_fs, data_fs)
    plt.savefig(model_dir + '_'.join(fs.pop_ids) + '_model_plot.png')
    plt.clf()

def make_demes_plot(model_dir, pop_ids, deme_model):
    '''
    Makes a demes plot from dadi model information
    Parameters:
        model_dir: String representing the model_specific dadi results directory
        pop_ids: List of strings of spp/pop names
        deme_model: Deme model output from Demo Model Optimization
    Returns:
        None
    '''
    deme_plot = demesdraw.tubes(deme_model)
    deme_plot.figure.savefig(model_dir + '_'.join(pop_ids) + '_demes_plot.png')
    plt.clf()

def godambe(popt, model_ex, pts, fs, model_dir, eps, boot_dir):
    '''
    Performs dadi Godambe Uncertainty Analysis on demographic models.
    Requires SFS bootstraps generated from dadi_2_sfs.py and will save the confidence intervals to text files in the dadi model results directory.
    Parameters:
        popt: List of int/floats of optimal model parameters
        model_ex: Our final model function
        pts: List of ints representing number of grid points for model
        fs: Data SFS object
        model_dir: String representing the model_specific dadi results directory
        eps: List of floats of step sizes to test confidence intervals
        boot_dir: String representing location of bootstrap data
    Returns:
        None
    '''
    # Get Bootstrapped datasets
    boots_fids = glob.glob(boot_dir + 'bootstraps/' + '_'.join(fs.pop_ids) + '/' + '_'.join(fs.pop_ids) + 'boots*.fs')
    boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

    # Godambe uncertainties will contain uncertainties for the estimated demographic parameters and theta.

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(fs.pop_ids) +'GIM_confidence_intervals.txt','w')
    fi.write('Optimized parameters: {0}\n\n'.format(popt))

    # We want to try a few different step sizes (eps) to see if uncertainties very wildly with changes to step size. (Ideally they shoud not)
    for steps in eps:
        uncerts_adj = dadi.Godambe.GIM_uncert(func_ex=model_ex, grid_pts=pts, all_boot=boots_syn, p0=popt, data=fs, eps=steps)
        uncerts_str = ',  '.join([str(ele) for ele in uncerts_adj])
        fi.write('GIM Array Output: [' + uncerts_str + ']\n')
        fi.write('Estimated 95% uncerts (with step size '+str(steps)+'): {0}\n'.format(1.96*uncerts_adj[:-1]))
        fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(popt-1.96*uncerts_adj[:-1]))
        fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(popt+1.96*uncerts_adj[:-1]))
            
    fi.close()
    
    

# Main
#==========================================================
'''
1) Check for required directories
2) Clean up CMD-Line Arguments
3) Generate demography model SFS and save fits to files
4) Generate Demes Plots
5) Plot SFS / Model SFS comparisons/residuals and save to files
6) Save uncertainty analysis on model to file
'''
#========================================
# Check if dadi-specific results directories exists in specified outdir. If not, create them.
print('Verifying Directories...')
# If using lowpass, set result dir to lowpass dir
result_dir = args.outdir + args.out_folder + '/lowpass/' if args.lowpass else args.outdir + args.out_folder + '/'
# The bootstrap files default to a nonlowpass directory, so we'll just store it in another variable to be used for GIM
boot_dir = args.outdir + args.out_folder + '/'
if not os.path.exists(result_dir):
    raise FileNotFoundError("Couldn't find specified result directory. Check that specified Result Directory matches the previously generated directory in data SFS creation.")

# Make new directory for model
model_dir = result_dir + args.model + '/'
if not os.path.exists(model_dir):
    os.makedirs(model_dir)

#========================================
# Clean up some arguments
up_params = [float(num) for num in json.loads(args.up_params[0])]
low_params = [float(num) for num in json.loads(args.low_params[0])]
start_params = [float(num) for num in json.loads(args.start_params[0])]
gim_steps = [float(step) for step in args.eps[0].split()]

#========================================
# Begin generating/optimizing model
print('Getting data fs for ' + args.job_name + '...')
fs = dadi.Spectrum.from_file(args.sfs_path)

# A Hack to determine fs dimensionality (If first element in fs array is NOT an array/spectrum object, then it is 1D)
# Thought there was a native dadi method to get dimensionality, but I guess not...
dim = 2 if type(fs[0]) == dadi.Spectrum_mod.Spectrum else 1

n = fs.sample_sizes
pts = [max(n)+20, max(n)+30, max(n)+40]

# State our model
if dim == 1:
    try:
        model = getattr(dadi.Demographics1D, args.model)
        print('Using 1d dadi Model')
    except:
        model = getattr(custom_models, args.model)
        print('Using 1d User Model')
else:
    try:
        model = getattr(dadi.Demographics2D, args.model)
        print('Using 2d dadi Model')
    except:
        model = getattr(custom_models, args.model)
        print('Using 2d User Model')

if not fs.folded:
    # Since we have unfolded data, we must wrap the model in a function that adds a parameter to estimate misidentification rate
    model = dadi.Numerics.make_anc_state_misid_func(model)

# Add another model wrapper that allows for the use of grid points
print('---> Make Model Function')
model_ex = dadi.Numerics.make_extrap_func(model)

if args.lowpass:
    from dadi.LowPass import LowPass
    import dill as pkl
    print('---> Loading Coverage Distribution...')
    with open(result_dir + '_'.join(fs.pop_ids) + '_cov_dist.pkl', 'rb') as file:
        cov_dist = pkl.load(file)
    print('---> Making LowPass Model Function...')
    model_ex = LowPass.make_low_pass_func_GATK_multisample(model_ex, cov_dist, fs.pop_ids, nseq=n, nsub=n, sim_threshold=1e-2, Fx=None)

# Create Fits File
fit_fname = model_dir + '_'.join(fs.pop_ids) + '_fits.txt'
with open(fit_fname,'w') as f:
    f.write('============ Optimized Model Fits ============\n')

# This is where we run the optimization for our arbitrary model parameters
# By the end, we will (presumably) get some optimal parameters and their log-likelihood
fits_lst = []

print('---> Perturb Starting Parameters...')
for i in range(args.num_opt):
    p0 = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=up_params,lower_bound=low_params)
    popt, ll_model = dadi.Inference.opt(p0, fs, model_ex, pts, lower_bound=low_params, upper_bound=up_params, algorithm=nlopt.LN_BOBYQA, maxeval=400, verbose=100, multinom=True)
    # Calculate the synonymous theta (Also finding optimal scaling factor for data)
    model_fs = model_ex(popt, n, pts)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

    # Write results to output file
    print('---> Saving Perturbed Model Parameters...')
    output = open(fit_fname,'a')
    res = [ll_model] + list(popt) + [theta0]
    fits_lst.append(res)
    output.write('\t'.join([str(ele) for ele in res])+'\n')
    output.close()

print('---> Finding Optimal Fits...')
# Prioritizing the parameter sets that have the best log-likelihood 
ll_opt = fits_lst[0][0]
fits_opt = fits_lst[0]
for lst in fits_lst:
    if lst[0] > ll_opt:
        ll_opt = lst[0]
        fits_opt = lst

print('---> Calculating Na...')
theta_opt = fits_opt[-1]
L = args.align_len * (fs.S() / args.num_snps)
Na = theta_opt / (4 * args.mu * L)
print(f"Na Calculation: Na = {theta_opt} / (4 * {args.mu} * ({args.align_len} * ({fs.S()} / {args.num_snps})))")
print('Na: ' + str(Na))

# If want to plot demes, do so here after model params are optimized
if args.demes:
    import demesdraw
    print('---> Saving Demes Plot...')
    # Rerun model with optimal params to properly store info for deme plotting
    model_opt = model_ex(fits_opt[1:-1], n, pts)
    deme_model_yrs = [dadi.Demes.output(Nref=Na, generation_time=args.generation_time), 'yrs']
    deme_model_gen = [dadi.Demes.output(Nref=Na), 'gen']
    deme_model_scale = [dadi.Demes.output(), 'scale']
    for d in [deme_model_yrs, deme_model_gen, deme_model_scale]:
        make_demes_plot(model_dir, fs.pop_ids, d)

print('Plotting SFS Comparison for ' + ' '.join(fs.pop_ids) + '...')
plot_model_fs(fs, model_fs, model_dir, dim)

print('Saving Model SFS file for ' + ' '.join(fs.pop_ids) + '...')
model_fs.to_file(model_dir + '_'.join(fs.pop_ids) + '_model_fs')

print('Running GIM Uncertainty Analysis')
godambe(popt, model_ex, pts, fs, model_dir, gim_steps, boot_dir)


print('\n**Model Creation Complete**')
