import dadi, nlopt, demesdraw
import matplotlib.pyplot as plt
import dill as pkl
import glob
import dfinbredmodels

dadi.cuda_enabled(True)

fs = dadi.Spectrum.from_file('dadi_results_cra_TEST_simple/CRA_pre_CRA_post_fs')

ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]
model = dfinbredmodels.df_split_mig
# SFS is folded, so we will not add misid param
model_ex = dadi.Numerics.make_extrap_func(model)

params = [0.01, 1, 1, 0.01, 0.5, 0.5]
lower_bounds = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]
upper_bounds = [10, 10, 10, 20, 1, 1]

try:
  fid = open('dadi_results_cra_TEST_simple/demo_fits.txt','a')
except:
  fid = open('dadi_results_cra_TEST_simple/demo_fits.txt','w')

for i in range(20):
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                              lower_bound=lower_bounds)
    popt, ll_model = dadi.Inference.opt(p0, fs, model_ex, pts,
                                    lower_bound=lower_bounds,
                                    upper_bound=upper_bounds,
                                    algorithm=nlopt.LN_BOBYQA,
                                    maxeval=400, verbose=100)
    model_fs = model_ex(popt, ns, pts)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)

res = [ll_model] + list(popt) + [theta0]
fid.write('\t'.join([str(ele) for ele in res])+'\n')
fid.close()

gim_params = [[popt, fs.pop_ids, model_ex, pts, fs]]

with open('dadi_results_cra_TEST_simple/gim_params.pkl', 'wb') as file:
        pkl.dump(gim_params, file)

fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(model_fs, fs, vmin=1)
fig.savefig('dadi_results_cra_TEST_simple/CRA_pre_CRA_post_model_plot.png')

deme_model = dadi.Demes.output(Nref=200)
deme_plot = demesdraw.tubes(deme_model)
deme_plot.figure.savefig('dadi_results_cra_TEST_simple/' + '_'.join(fs.pop_ids) + '_demes_plot.png')
plt.clf()

def load_bootstraps(result_dir, pop_ids):
    # Get Bootstrapped datasets
    boots_fids = glob.glob(result_dir + 'bootstraps/' + '_'.join(pop_ids) + '/' + '_'.join(pop_ids) + 'boots*.fs')
    boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]
    return boots_syn

def write_to_file(fi, uncerts, steps, popt):
    uncerts_str = ',  '.join([str(ele) for ele in uncerts])
    fi.write('GIM Array Output: [' + uncerts_str + ']\n')
    fi.write('Estimated 95% uncerts (with step size '+str(steps)+'): {0}\n'.format(1.96*uncerts[:-1]))
    fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(popt-1.96*uncerts[:-1]))
    fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(popt+1.96*uncerts[:-1]))

def godambe(popt, model_ex, pts, fs, model_dir, eps, result_dir, pop_ids):
    # Get Bootstrapped datasets
    boots_syn = load_bootstraps(result_dir, pop_ids)

    # Godambe uncertainties will contain uncertainties for the estimated demographic parameters and theta.

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(pop_ids) +'GIM_confidence_intervals.txt','w')
    fi.write('Optimized parameters: {0}\n\n'.format(popt))

    # We want to try a few different step sizes (eps) to see if uncertainties very wildly with changes to step size. (Ideally they shoud not)
    for steps in eps:
        uncerts_adj = dadi.Godambe.GIM_uncert(func_ex=model_ex, grid_pts=pts, all_boot=boots_syn, p0=popt, data=fs, eps=steps, log=True)
        write_to_file(fi, uncerts_adj, steps, popt)
    fi.close()

def likelihood(popt, model_ex, pts, fs, model_dir, eps, result_dir, pop_ids, lrt_indices):
    # Get Bootstrapped datasets
    boots_syn = load_bootstraps(result_dir, pop_ids)

    # Start a file to contain the confidence intervals
    fi = open(model_dir  + '_'.join(pop_ids) +'LRT_results.txt','w')
    fi.write('Optimized parameters: {0}\n\n'.format(popt))

    for steps in eps:
        adj = dadi.Godambe.LRT_adjust(func_ex=model_ex, grid_pts=pts, all_boot=boots_syn, p0=popt, data=fs, nested_indices=lrt_indices, multinom = True, eps=steps)
        write_to_file(fi, adj, steps, popt)
    fi.close()
    
godambe(popt, model_ex, pts, fs, 'dadi_results_cra_TEST_simple/', [0.01,0.001, 0.0001, 0.00001], 'dadi_results_cra_TEST_simple/', fs.pop_ids)

#likelihood(popt, model_ex, pts, fs, 'dadi_results_cra_TEST_simple/', [0.01,0.001, 0.0001, 0.00001], 'dadi_results_cra_TEST_simple/', fs.pop_ids, [3])
