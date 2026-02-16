# LDHelmet Parameters

micromamba activate ldhelmet_env # If set up using a conda/mamba environment
source ../../../params_base.sh

INPUT_DIR=${OUTDIR}/datafiles/ldhelmet_inputs/
RESULT_DIR=${OUTDIR}/analyses/ldhelmet_outputs/

# General Parameters
THREADS=<SET_VALUE>
WINDOW_SIZE=<SET_VALUE> # often 50
MUT_RATE=<SET_VALUE> # Population scaled mutation rate in units of 1/bp

# Likelihood Lookup Table
REC_RATE_GRID=0.0 0.1 10.0 1.0 100.0

# Pade
PADE_COEF=<SET_VALUE>
DEFECT=<SET_VALUE>

# MCMC Analysis
BURN_IN=<SET_VALUE>
BLOCK_PENALTY=<SET_VALUE>
ITERATIONS=<SET_VALUE>
MUT_MATRIX= # Optional, leave empty if not using a custom mutation matrix. If using, provide the path to the mutation matrix file.
ANC_PRIOR= # Optional, leave empty if not using ancestral allele priors. If using, provide the path to the ancestral allele prior file. 