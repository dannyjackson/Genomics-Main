'''
Author: <Logan Vossler>

==========================================================
Run Likelihood-Ratio Test
==========================================================

Description:
Description:
This script uses dadi to perform a Likelihood Ratio Test between two models.

Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - A Data SFS File
    - Bootstrapped SFS files in proper result directory
    - The two models that you are comparing must be "nested". Meaning that they should be very similar models that differ in only some parameters;
        For example: You could compare a simple bottlegrowth model to a more complex bottlegrowth model that has an additional inbreeding param.
        See dadi docs for additional info on requirements for models to be nested.
'''


# Required Modules
#==========================================================
import dadi, glob, os
import dfinbredmodels as custom_models
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-j", "--job_name", type=str, help='Name of Job')
parser.add_argument("-n", "--null_model", type=str, help='Null Demographic Model')
parser.add_argument("-s", "--sfs_path", type=str, help='Path to Data SFS')
parser.add_argument("-o","--outdir", type=str, help='Out Directory Path')
parser.add_argument("-t","--test_model", type=str, help='Test Demographic Model')
parser.add_argument("-i", "--nested_idx", type=str, nargs='+', help='List of indices')
parser.add_argument("--null_popt", type=str, nargs='+', help='List of optimized null model params')
parser.add_argument("--test_popt", type=str, nargs='+', help='List of optimized test model params')
parser.add_argument("-b", "--boot_dir", type=str, help='Path to Bootstrap Directory')
args = parser.parse_args()

# Main
#==========================================================
'''
1) Import Base Parameters
2) Import Dadi-SFS Specific Parameters
3) Check for required directories
4) Load GIM Parameters
5) Perform GIM Uncertainty Analysis and write results to files
'''
#========================================
# Clean up some arguments
test_model = args.test_model.split('.')[-1]
null_model = args.null_model.split('.')[-1]
nest_idx = [int(idx) for idx in args.nested_idx[0].split()]
null_popt = [float(param) for param in args.null_popt[0].split()]
test_popt = [float(param) for param in args.test_popt[0].split()]


#========================================
# Check if dadi-specific results directories exists in specified outdir. If not, create them.
print('Verifying Directories...')
result_dir = args.outdir + 'LRT_outputs/'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

#========================================
print('Loading SFS...')
data_fs = dadi.Spectrum.from_file(args.sfs_path)
n = data_fs.sample_sizes
pts = [max(n)+20, max(n)+30, max(n)+40]


print('Building Models...')
# Our model dimensionality should be the same as our data dimensionality, therefore...
model_dim = 2 if type(data_fs[0]) == dadi.Spectrum_mod.Spectrum else 1

if model_dim == 2:
    try:
        test_func = getattr(dadi.Demographics2D, test_model)
    except:
        test_func = getattr(custom_models, test_model)

    try:
        null_func = getattr(dadi.Demographics2D, null_model)
    except:
        null_func = getattr(custom_models, null_model)

elif model_dim == 1:
    try:
        test_func = getattr(dadi.Demographics1D, test_model)
    except:
        test_func = getattr(custom_models, test_model)

    try:
        null_func = getattr(dadi.Demographics1D, null_model)
    except:
        null_func = getattr(custom_models, null_model)

else:
    print(f"Model Dimensionality: {model_dim}")
    raise ValueError('This script only supports either 1 or 2 dimensional SFS and models!')


if not data_fs.folded:
    test_func = dadi.Numerics.make_anc_state_misid_func(test_func)
    null_func = dadi.Numerics.make_anc_state_misid_func(null_func)


test_ex = dadi.Numerics.make_extrap_log_func(test_func)
null_ex = dadi.Numerics.make_extrap_log_func(null_func)



print('Finding Model Log-Likelihoods...')
model_test = test_ex(test_popt, n, pts)
model_null = null_ex(null_popt, n, pts)

ll_test = dadi.Inference.ll_multinom(model_test, data_fs)
ll_null = dadi.Inference.ll_multinom(model_null, data_fs)

print(f"Test Model ll: {ll_test} \nNull Model ll: {ll_null}")


#=============================================================
# Time for LRT

# Since LRT evaluates the complex (test) model using the best-fit parameters from the
# simple (null) model, we need to create a list of parameters for the complex model
# using the simple best-fit params.  Since evalution is done with the more
# complex model, we need to insert zero values at corresponding
# parameter indices in complex model. And we need to tell the LRT adjust function
# that these parameters (counting from 0) are the nested ones.

# Grab Bootstraps
boots_fids = glob.glob(os.path.join(args.boot_dir, '*'))
boot_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

# Insert zero-ed params to null_popt lst
for idx in nest_idx:
    null_popt.insert(idx, 0)

adj = dadi.Godambe.LRT_adjust(test_ex, pts, boot_syn, null_popt, data_fs, nested_indices=nest_idx, multinom=True)
print('adj: ' + str(adj))

D_adj = adj*2*(ll_test - ll_null)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))
Dstat_str = 'Adjusted D statistic: {0:.4f}'.format(D_adj)

# Because this is a test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting null model: {0:.4f}'.format(pval))
pval_str = 'p-value for rejecting null model: {0:.4f}'.format(pval)

with open(os.path.join(result_dir, test_model + '_' + null_model + 'LRT.out'), 'w') as file:
    file.write(Dstat_str)
    file.write(pval_str)
