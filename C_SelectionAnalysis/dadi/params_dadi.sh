# dadi Specific Parameters

# Universal Params
# Parameters used generally across all or multiple analyses
JOB_NAME=par_post_1d_model_generation
OUT_FOLDER=dadi_results_par # The name of the folder (under your general OUTDIR) you want to generate containing your results NOTE: Don't forgot to update this when switching the populations being analyzed (otherwise the scripts will not output your results in desired location or find needed intermediate param files)
LOWPASS=TRUE # Must be set to True during SFS creation so that the lowpass coverage distribution is generated (Though this is awkward, we generate the cov-dist during SFS creation so that way our model scripts don't need to waste resources reloading data dictionaries)
SFS_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/dadi_results_par/PAR_post_fs' # Used in model generations and running LRT analysis


# SFS Creation Params
VCF_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/vcfs/cra_all_qualitysort.vcf'
POP_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/vcfs/cra_post_pops.txt'
POP_IDS="PAR_post" # Space Separated populations IDs being analyzed
NUM_CHROMS="18" # Space separated numbers representing number of chromosomes per population
BOOTSTRAP_PARAMS="100 1e-7" # Space separated numbers representing the number of bootstraps to be generated and the chunk size respectively
POLARIZE=FALSE # Switch this to TRUE if you want an UNFOLDED (not triangular) SFS (meaning you are confident of the identities of the ancestral/derived alleles)


# Model Creation Params
NUM_OPT=100
PLOT_DEMES=TRUE
GIM_STEPS="0.1 0.01 0.001 0.0001 0.00001" # Space Separated list of step sizes
MODEL_JSON='/xdisk/mcnew/finches/ljvossler/finches/dadi/param_files/dadi_1dmodel_params.json'


# LRT Analysis Params
NESTED_INDICES="4 5" #Should be a space separated list of indices 
TEST_MODEL="dfinbredmodels.split_mig_inbred"
NULL_MODEL="dadi.Demographics2D.split_mig"
NULL_POPT="1.52155171641276 0.08679806754804995 0.001 0.005963741213696863"    # Space separated list of optimized model parameters (NOTE: You must include zero-ed param valuesin indices to be used in complex model
TEST_POPT="3.0 3.0 0.001 0.02118775400887089 0.001 0.23288339760430127"    # Space separated list of optimized model parameters
BOOT_DIR='/xdisk/mcnew/finches/ljvossler/finches/dadi/dadi_results_cra/bootstraps/CRA_pre_CRA_post/'