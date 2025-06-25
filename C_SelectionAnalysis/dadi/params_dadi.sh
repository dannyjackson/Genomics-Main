# dadi Specific Parameters

# Universal Params
# Parameters used generally across all analyses
JOB_NAME=par_sfs_generation
OUT_FOLDER=dadi_results_unfiltered_par # The name of the folder (under your general OUTDIR) you want to generate containing your results NOTE: Don't forgot to update this when switching the populations being analyzed (otherwise the scripts will not output your results in desired location or find needed intermediate param files)
POP_IDS="PAR_pre PAR_post" # Space Separated populations IDs being analyzed
LOWPASS=TRUE # Must be set to True during SFS creation so that the lowpass coverage distribution is generated (Though a bit awkward, we generate the cov-dist during SFS creation so that way our model scripts don't need to waste resources reloading data dictionaries)


# SFS Creation Params
VCF_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/vcfs/par_all_qualitysort.vcf'
POP_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/vcfs/par_pops.txt'
NUM_CHROMS="24 14" # Space separated numbers representing number of chromosomes per population
BOOTSTRAP_PARAMS="100 1e-7" # Space separated numbers representing the number of bootstraps to be generated and the chunk size respectively
POLARIZE=FALSE # Switch this to TRUE if you want an UNFOLDED (not triangular) SFS (meaning you are confident of the identities of the ancestral/derived alleles)


# Model Creation Params
NUM_OPT=20
PLOT_DEMES=TRUE
SFS_PATH='/xdisk/mcnew/finches/ljvossler/finches/dadi/dadi_results_unfiltered_cra/CRA_pre_CRA_post_fs'
GIM_STEPS="0.1 0.01 0.001 0.0001 0.00001"
MODEL_JSON='/xdisk/mcnew/finches/ljvossler/finches/dadi/param_files/test_params.json'

# LRT Analysis Params
NESTED_INDICES=3 #Should be a space separated list of indices 
SFS_PATH=/path/to/data/sfs
TEST_MODEL="dadi.Demographics2D.split_mig"
NULL_MODEL="dadi.Demographics2D.snm_2d"
NULL_POPT=             # Space separated list of optimized model parameters
TEST_POPT=              # Space separated list of optimized model parameters
BOOTSTRAP_DIR=/path/to/bootstraps/dir