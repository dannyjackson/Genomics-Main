# flexsweep_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/analyses/flexsweep"

echo "Installing flexsweep executable with uv to user home directory"
uv tool install flexsweep

echo "FlexSweep setup completed."