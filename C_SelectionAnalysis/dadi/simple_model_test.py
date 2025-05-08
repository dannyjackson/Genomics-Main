import dadi, nlopt
import matplotlib.pyplot as plt

dadi.cuda_enabled(True)

fs = dadi.Spectrum.from_file('dadi_results_unfiltered_cra/lowpass/CRA_pre_CRA_post_fs')

ns = fs.sample_sizes
pts = [max(ns)+20, max(ns)+30, max(ns)+40]
model = dadi.Demographics2D.split_mig
# SFS is folded, so we will not add misid param
model_ex = dadi.Numerics.make_extrap_func(model)

params = [1, 1, 0.01, 0.01, 0.01]
lower_bounds = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3]
upper_bounds = [3, 3, 1, 1, 1]

try:
  fid = open('cuda_test_demo_fits.txt','a')
except:
  fid = open('cuda_test_demo_fits.txt','w')

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

fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(model_fs, fs)
fig.savefig('cuda_test_demo_plot.png')