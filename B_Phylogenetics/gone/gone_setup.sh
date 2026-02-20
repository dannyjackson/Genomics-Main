# gone2_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/datafiles/gone2_inputs"
mkdir -p "${OUTDIR}/analyses/gone2_outputs"

echo "Getting GONE2 scripts from Github..."
date

# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d GONE2 ]; then
    git clone https://github.com/esrud/GONE2.git || { echo "Error: Failed to clone GONE2."; exit 1; }
fi

if [ -d "GONE2" ]; then
    echo "GONE2 directory exists."
    cd GONE2
    make gone # Edit this as needed based on your system and project memory requirements. See GONE2 documentation for details.
else
    echo "Error: Failed to build GONE2."
fi

echo "GONE2 Environment setup completed."