# LDhelmet_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/datafiles/ldhelmet"

echo "Getting LDHelmet source from Github..."
date

# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d GONE2 ]; then
    git clone https://github.com/popgenmethods/LDhelmet.git || { echo "Error: Failed to clone GONE2."; exit 1; }
fi

# Getting LDHelmet examples from github is useful for testing with the authors' examples. But on the UA HPC, it is easiest to just install it using micromamba.

#LDHelmet is on Bioconda, and its required dependencies (Boost and GSL) are on conda-forge, so ensure your micromamba config includes these required channels.
module load micromamba
micromamba config append channels conda-forge
micromamba config append channels bioconda

micromamba config list

# Install LDHelmet using micromamba in a new environment
micromamba create -n ldhelmet_env ldhelmet


echo "LDHelmet Environment setup completed."