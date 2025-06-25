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
import dfinbredmodels

print('define arguments')
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

# Function Definitions
#==========================================================
def load_bootstraps(boot_dir):
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
    test_model = args.test_model.split('.')[-1]
    null_model = args.null_model.split('.')[-1]
    nest_idx = args.nested_idx[0].split()
    null_popt = args.null_popt[0].split()
    test_popt = args.test_popt[0].split()
    

    #========================================
    # Check if dadi-specific results directories exists in specified outdir. If not, create them.
    print('Verifying Directories...')
    result_dir = args.outdir + 'LRT_outputs/'
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

    # Get data_fs
    data_fs = dadi.Spectrum.from_file(args.sfs_path)

    # Get log-likelihoods for each model
    n = data_fs.sample_sizes
    pts = [max(n)+20, max(n)+30, max(n)+40]
    model_test = test_ex(test_popt, n, pts)
    model_null = null_ex(null_popt, n, pts)
    ll_test = dadi.Inference.ll_multinom(model_test, data_fs)
    ll_null = dadi.Inference.ll_multinom(model_null, data_fs)

    # Since LRT evaluates the complex model using the best-fit parameters from the
    # simple model, we need to create list of parameters for the complex model
    # using the simple (no-mig) best-fit params.  Since evalution is done with more
    # complex model, need to insert zero migration value at corresponding migration
    # parameter index in complex model. And we need to tell the LRT adjust function
    # that the 3rd parameter (counting from 0) is the nested one.

    # Grab Bootstraps
    boot_syn = load_bootstraps(args.boot_dir)

    adj = dadi.Godambe.LRT_adjust(test_ex, pts, boot_syn, null_popt, data_fs, nested_indices=nest_idx, multinom=True)
    D_adj = adj*2*(ll_test - ll_null)
    print('Adjusted D statistic: {0:.4f}'.format(D_adj))
    Dstat_str = 'Adjusted D statistic: {0:.4f}'.format(D_adj)

    # Because this is test of a parameter on the boundary of parameter space 
    # (m cannot be less than zero), our null distribution is an even proportion 
    # of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
    # point percent function for a weighted sum of chi^2 dists.
    pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
    print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))
    pval_str = 'p-value for rejecting no-migration model: {0:.4f}'.format(pval)

    with open(os.path.join(result_dir, test_model + '_' + null_model + 'LRT.out'), 'w') as file:
        file.write(Dstat_str)
        file.write(pval_str)
