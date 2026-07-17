# preprocessing_setup.sh

# A0

# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/condensed_fastas" "${OUTDIR}/datafiles/trimming"
mkdir -p "${OUTDIR}/datafiles/trimmed_fastas" "${OUTDIR}/datafiles/bamfiles/${IND}" "${OUTDIR}/datafiles/sortedbamfiles/${IND}"
mkdir -p "${OUTDIR}/datafiles/clipoverlap"
mkdir -p "${OUTDIR}/datafiles/indelmaps"
mkdir -p "${OUTDIR}/datafiles/bamstats"
mkdir -p "${OUTDIR}/datafiles/indelrealignment"

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

if [ ! -d sapphire ]; then
    git clone https://github.com/stschiff/sapphire || { echo "Error: Failed to clone sapphire."; exit 1; }
fi
cd sapphire # Note that you should be on the Puma Cluster to successfully build and run sapphire
module load htslib
make
cd ../
#======================================================

if [ ! -d pyrho ]; then
    echo "Installing pyrho..." # Note that you should be on the Puma Cluster to successfully build and run pyrho

    micromamba create -y -n pyrho_env \
            python=3.7 \
            numpy=1.19.* \
            scipy=1.5.* \
            pandas=1.1.* \
            msprime=1.0.* \
            numba=0.53.* \
            pytables=3.6.* \
            cyvcf2=0.11.* \
            htslib=1.9.* \
            pytest \
            -c conda-forge -c bioconda

    git clone https://github.com/popgenmethods/ldpop.git ldpop
    pip install ldpop/
    git clone https://github.com/popgenmethods/pyrho.git pyrho
    pip install pyrho/

    python -m pytest ${PROGDIR}/pyrho/tests/tests.py # Ensure that all tests are passed
fi
# So there is an issue with the current commit of pyrho (as of July 2026). This may just be the way we are installing on the HPC, not necessarily an issue with pyrho itself.
# Running hyperparams throws an error since a datafile argument in the code is invalid. 
# This can be fixed my manually editing the code where it throws the error (remove the subdirectory path), reinstall pyrho, then move that file in the environments folder up one directory
# ****come back to this later to see if there is just a minor issue in installing that I should fix****

#======================================================

echo "Installing SMC++..."
apptainer pull docker://terhorst/smcpp
#======================================================

echo "Installing BEAGLE5.5..."
micromamba create -n beagle_env -c bioconda beagle
#======================================================