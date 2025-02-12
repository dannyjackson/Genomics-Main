# load required modules
module load gnu13/13.2.0
module load gsl/2.7.1
module load R/4.4.0

# make directories for intermediate files-- will fail if these don't exist
mkdir -p ${OUTDIR}/analyses/raisd/
mkdir -p ${OUTDIR}/analyses/raisd/${POP}
mkdir -p ${OUTDIR}/analyses/raisd/${POP}/${WIN}
mkdir -p ${OUTDIR}/analyses/raisd/${POP}/${WIN}/{infofiles,plots,reportfiles}


# install raisd
RAiSD_DIR="${PROGDIR}/RAiSD"
RAiSD_EXEC="${RAiSD_DIR}/raisd-master/RAiSD"

# Check if RAiSD is already installed
if [ -f "$RAiSD_EXEC" ]; then
    echo "RAiSD is already installed."
else
    echo "Installing RAiSD..."
    
    # Create directory and download RAiSD
    mkdir -p "$RAiSD_DIR"
    cd "$RAiSD_DIR" || exit
    wget https://github.com/alachins/raisd/archive/master.zip
    unzip master.zip
    cd raisd-master || exit
    
    # Install RAiSD
    ./install-RAiSD.sh
    
    # Update LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${RAiSD_DIR}/raisd-master/gsl/lib
    
    echo "RAiSD installation complete."
fi