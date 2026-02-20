# base_setup.sh

# make main directories
# specific to selection analyses (fst, dxy, Tajima's D, RAiSD)
# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/msmc/outputs/
mkdir -p ${OUTDIR}/analyses/msmc/bootstraps/
mkdir -p ${OUTDIR}/datafiles/msmc/input/
mkdir -p ${OUTDIR}/datafiles/msmc/mask/ind/
mkdir -p ${OUTDIR}/datafiles/msmc/mask/genom/
mkdir -p ${OUTDIR}/datafiles/bootstraps/boots_file_lsts/


# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d msmc2 ]; then
    git clone https://github.com/stschiff/msmc2 || { echo "Error: Failed to clone msmc2."; exit 1; }
fi

if [ ! -d msmc-tools ]; then
    git clone https://github.com/stschiff/msmc-tools || { echo "Error: Failed to clone msmc-tools."; exit 1; }
fi

# The UA HPC does not seem to support a D Compiler, which is required for building MSMC2 from source. So we will just use precompiled exectutable from bioconda (https://bioconda.github.io/recipes/msmc2/README.html)

module load micromamba
micromamba create -n msmc_env python=3.11 # Add python to environment for msmc-tools
micromamba micromamba activate msmc_env
micromamba install msmc2

