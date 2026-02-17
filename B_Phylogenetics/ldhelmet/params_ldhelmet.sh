# LDHelmet Parameters

module load micromamba
micromamba activate ldhelmet_env # If set up using a conda/mamba environment
source ../../../params_base.sh

INPUT_DIR=${OUTDIR}/datafiles/ldhelmet

# General Parameters
THREADS=12
WINDOW_SIZE=50 # often 50
MUT_RATE=<SET_VALUE> # Population scaled mutation rate in units of 1/bp

# Likelihood Lookup Table
REC_RATE_GRID=0.0 0.1 10.0 1.0 100.0

# Pade
PADE_COEF=12
DEFECT=40

# MCMC Analysis
BURN_IN=10000
BLOCK_PENALTY=50
ITERATIONS=100000
MUT_MATRIX= # Optional, leave empty if not using a custom mutation matrix. If using, provide the path to the mutation matrix file.
ANC_PRIOR= # Optional, leave empty if not using ancestral allele priors. If using, provide the path to the ancestral allele prior file. 