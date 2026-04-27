# preprocessing_setup.sh

# A0

# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/condensed_fastas" "${OUTDIR}/datafiles/trimming"
mkdir -p "${OUTDIR}/datafiles/trimmed_fastas" "${OUTDIR}/datafiles/bamfiles/${IND}" "${OUTDIR}/datafiles/sortedbamfiles/${IND}"
mkdir -p "${OUTDIR}/datafiles/clipoverlap"
mkdir -p "${OUTDIR}/datafiles/indelmaps"
mkdir -p "${OUTDIR}/datafiles/bamstats"
mkdir -p "${OUTDIR}/datafiles/indelrealignment"
mkdir -p "${OUTDIR}/datafiles/sapphire_polishing"

# A1
# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/geno_likelihoods"
mkdir -p "${OUTDIR}/datafiles/genotype_calls/"
mkdir -p "${OUTDIR}/datafiles/snpable"
mkdir -p "${OUTDIR}/datafiles/vcf"
mkdir -p "${OUTDIR}/datafiles/vcf2"
mkdir -p "${OUTDIR}/datafiles/mask"
mkdir -p "${OUTDIR}/datafiles/stats"
mkdir -p "${OUTDIR}/datafiles/safs"

# A2 
# Define path for Snpable scripts
SNPABLE_SCRIPT_PATH="${PROGDIR}/seqbility-20091110"

# Check if the Snpable directory exists
if [ ! -d "$SNPABLE_SCRIPT_PATH" ]; then
    echo "Error: Snpable script directory not found at ${SNPABLE_SCRIPT_PATH}." >&2
    exit 1
fi
# Add Snpable scripts to PATH
export PATH="$PATH:$SNPABLE_SCRIPT_PATH"



echo "Installing Various Program dependencies..."


#======================================================
echo "Adding tools for MSMC pipeline..."
date
# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d msmc2 ]; then
    git clone https://github.com/stschiff/msmc2 || { echo "Error: Failed to clone msmc2."; exit 1; }
fi
if [ ! -d msmc-tools ]; then
    git clone https://github.com/stschiff/msmc-tools || { echo "Error: Failed to clone msmc-tools."; exit 1; }
fi
#======================================================

echo "Installing whatshap..."
pip3 install --user whatshap
#======================================================

echo "Installing shapeit4..."
micromamba create -n shapeit4_env -c bioconda shapeit4

if [ ! -d sapphire ]; then
    git clone https://github.com/stschiff/sapphire || { echo "Error: Failed to clone sapphire."; exit 1; }
fi
cd sapphire # Note that you should be on the Puma Cluster to successfully build and run sapphire
module load htslib
make
cd ../
#======================================================

echo "Installing pyrho..." # Note that you should be on the Puma Cluster to successfully build and run pyrho
micromamba create -n pyrho_env
micromamba activate pyrho_env
micromamba install python=3.11 numpy=2.4 cython
module load gsl
module load htslib
module load hdf5
git clone https://github.com/popgenmethods/ldpop.git ldpop
pip install ldpop/
pip install cython
git clone https://github.com/popgenmethods/pyrho.git pyrho
pip install pyrho/
#======================================================

echo "Installing SMC++..."
apptainer pull docker://terhorst/smcpp
#======================================================

echo "Installing BEAGLE5.5..."
micromamba create -n beagle_env -c bioconda beagle
#======================================================