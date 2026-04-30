# flexsweep_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/analyses/flexsweep_outputs"

echo "Getting FlexSweep using micromamba..."
date

micromamba create -n flexsweep_env -c bioconda flexsweep # ensure you are installing while on Puma. You'll need to run the program from here as well.

echo "FlexSweep Environment setup completed."