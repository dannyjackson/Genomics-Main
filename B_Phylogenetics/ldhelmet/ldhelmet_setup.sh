# LDhelmet_setup.sh
# Create necessary directories

source ../../../params_base.sh

mkdir -p "${OUTDIR}/datafiles/ldhelmet"
mkdir -p "${OUTDIR}/analyses/ldhelmet"

echo "Getting LDHelmet source from Github..."
date

#======================================
#IMPORTANT!!!!

# It is HIGHLY recommended to install and compile LDHelmet from its git repo. 

# You may find it easier to install LDHelmet from bioconda using micromamba. However, you should be aware that the final post-processing steps for LDHelmet are very memory intensive and LDHelmet has a HARD-CODED memory limit of 20gb 
# This cannot be changed with any simple component flags. In this scenario you'll have to be prepared to split your data into smaller partitions to reduce memory load.

# Therefore, to avoid this for an easier time in the long-run I'd strongly recommend installing LDHelmet from github and making some required edits to the source code and Makefile.
#======================================


# Clone required repositories if not already present
cd "${PROGDIR}" || { echo "Error: Could not change directory to ${PROGDIR}."; exit 1; }

if [ ! -d GONE2 ]; then
    git clone https://github.com/popgenmethods/LDhelmet.git || { echo "Error: Failed to clone GONE2."; exit 1; }
fi

echo "To build LDHelmet, refer to the instructions in Sections 2.3 and 2.4 of the LDHelmet manual"

printf "\n\n"

echo "Your *include* and *lib* paths for GSL are under the prepend_paths below:"
module show gsl

printf "\n\n"

echo "Your *include* and *lib* paths for Boost are under the prepend_paths below:"
module show boost

printf "\n\n"

#=======================================
echo "To remove the memory limit in LDHelmet, open 'src/common/site_map_log_lk.cc' in the LDHelmet git directory and remove/comment out the lines that read:
    
    if (size_site_log_lks > 2684354560) {
    fprintf(stderr,
            "Safety check: "
            "You probably don't want to analyze such a large number "
            "of SNPs at once. The amount of memory required will be "
            "at least 20 GB. If you do, you'll have to remove this "
            "safety check from the code.\n"
            "Recommendation: Use shorter partitions.\n");
    std::exit(1);
  }
    "
#=======================================

# Alternative micromamba installation
#module load micromamba
#micromamba config append channels conda-forge
#micromamba config append channels bioconda
#micromamba config list
#micromamba create -n ldhelmet_env ldhelmet
