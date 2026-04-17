#!/bin/bash


# THIS SCRIPT IS UNTESTED!!!!!!!!!!!!!!

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-s]"
    echo ""
    echo "This script phases a VCF file with BEAGLE5.5 using a VCF file containing all individuals in a desired population. It can also polish the phased output with SAPPHIRE, which is highly recommended"
    echo "It is best run as a Slurm array that calls this script for each individual/population."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the population to analyze."
    echo "  -s  Optional flag to enable SAPPHIRE polishing of phased VCFs (recommended)."
    exit 1
}

SAPPHIRE=false

# Parse command-line arguments
while getopts ":p:i:s" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        s) SAPPHIRE=true ;; # False by default, recommended to set to true
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$IND" || -z "$BAMDIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Load parameters
source "${PARAMS}"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

VCF_IN="${OUTDIR}/datafiles/vcf/${POP}_all.vcf.gz"
VCF_OUT="${OUTDIR}/datafiles/phased_vcf/${POP}_phased.vcf.gz"

beagle gt=${VCF_IN} out=${VCF_OUT} map=${OUTDIR}/datafiles/linkage_map/${POP}_all.rmap

bcftools index ${VCF_OUT}

# Perform sapphire rephasing

if [ "$SAPPHIRE" = true ]; then
    ${SCRIPT_DIR}/Genomics-Main/A_Preprocessing/A2.6_rephasing.sh -p "${PARAMS}" -i "${POP}" -v "${VCF_OUT}"


echo "VCF phasing completed."
date
