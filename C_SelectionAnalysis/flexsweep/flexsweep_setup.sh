# flexsweep_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/analyses/flexsweep"

#echo "Getting FlexSweep using micromamba..."
#micromamba create -n flexsweep_env -c bioconda flexsweep

echo "Installing flexsweep with uv"
cd ${PROGDIR}
uv flexsweep_venv --python 3.12
source .flexsweep_venv/bin/activate
uv pip install flexsweep

echo "FlexSweep Environment setup completed. Environment located in ${PROGDIR}/flexsweep_venv/"