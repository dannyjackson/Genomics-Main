#!/bin/sh

# raisd script

# RAiSD script for detecting selective sweeps from whole-genome sequence data.

# Function to display usage information
usage() {
    cat <<EOF
This script applies the software RAiSD to detect selective sweeps from whole-genome sequence data.

It requires a gzipped vcf file as input, which can be generated using the <scriptname> scripts in:
    github.com/dannyjackson/Genomics-Main

REQUIRED ARGUMENTS:
    -w  Window size
    -p population name
EOF
    exit 1
}


# Parse command-line arguments
while getopts "p:w:" option; do
    case "${option}" in
        p) POP=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done


# Source the parameter file 
source ~/programs/CardinalisGenomics/params_base.sh

# Ensure that necessary parameters are set in the parameter file
if [ -z "${OUTDIR:-}" ] ]; then
    echo "Error: OUTDIR is not set in the parameter file." >&2
    exit 1
fi


# Print starting info
echo -e "\n\n\n\n"
date
echo "Current script: raisd.sh"

# Ensure the output directory exists and navigate there
cd "${OUTDIR}/analyses/raisd/${POP}" || { echo "Failed to enter directory"; exit 1; }

# Validate that the scaffold list file exists
SCAFFOLD_LIST="${OUTDIR}/referencelists/SCAFFOLDS.txt"
if [[ ! -f "$SCAFFOLD_LIST" ]]; then
    echo "Error: Scaffold list file not found: $SCAFFOLD_LIST"
    exit 1
fi

# Process each scaffold
while read -r SCAFFOLD; do
    echo "Processing scaffold: $SCAFFOLD"
    
    # Define input VCF file path
    VCF_IN="${OUTDIR}/datafiles/vcf2/${POP}.${SCAFFOLD}.phased.vcf.gz"
    
    # Check if the VCF file exists
    if [[ ! -f "$VCF_IN" ]]; then
        echo "Error: VCF file not found: $VCF_IN"
        continue
    fi
    
    # Run RAiSD program (ensure the path to RAiSD is correct)
    ${PROGDIR}/RAiSD/raisd-master/RAiSD \
        -n "${OUTDIR}/raisd/output/${POP}.${SCAFFOLD}" \
        -I "${VCF_IN}" \
        -f -O -R -P -a 1500 -C "${REF}" -w ${WIN}

done < "$SCAFFOLD_LIST"

# Completion message
echo "RAiSD analysis completed for all scaffolds."
date