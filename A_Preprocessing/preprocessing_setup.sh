# preprocessing_setup.sh

# A0

# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/condensed_fastas" "${OUTDIR}/datafiles/trimming"
mkdir -p "${OUTDIR}/datafiles/trimmed_fastas" "${OUTDIR}/datafiles/bamfiles/${IND}" "${OUTDIR}/datafiles/sortedbamfiles/${IND}"
mkdir -p "${OUTDIR}/datafiles/clipoverlap"
mkdir -p "${OUTDIR}/datafiles/indelmaps"
mkdir -p "${OUTDIR}/datafiles/bamstats"

# A1
# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/geno_likelihoods"
mkdir -p "${OUTDIR}/datafiles/genotype_calls/"
mkdir -p "${OUTDIR}/snpable"

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

echo "Starting MSMC pipeline..."
date

# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d msmc2 ]; then
    git clone https://github.com/stschiff/msmc2 || { echo "Error: Failed to clone msmc2."; exit 1; }
fi

if [ ! -d msmc-tools ]; then
    git clone https://github.com/stschiff/msmc-tools || { echo "Error: Failed to clone msmc-tools."; exit 1; }
fi

echo "Environment setup completed."


# Install required dependencies
pip3 install --user whatshap
