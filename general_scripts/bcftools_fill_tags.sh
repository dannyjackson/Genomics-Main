#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script filters the VCF generated in A2.1_callvariants.sh"
    echo
    echo "Required Arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh from the GitHub repository)"
    echo "  -v  Path to VCF file"
    exit 1
}

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    usage
fi

# Parse command-line arguments
while getopts "p:v:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        v) VCF=${OPTARG} ;;
        *) usage ;;
    esac
done

source ${PARAMS}

echo Inputted VCF: $VCF
BCF="${VCF%.vcf.gz}".bcf

bcftools view $VCF -Ob -o ${BCF}
bcftools index ${BCF}

# Fill extra tags in filtered VCF. Currently this command will populate allele frequency, count, and number fields. 
# We are also converting to BCFs here. Not required, but good for SAPPHIRE rephasing (A2.7_rephasing.sh).
bcftools +fill-tags $BCF -Ob -o "${BCF%.bcf}"_tagfilled.bcf -- -t AF,AC,AN
bcftools index "${BCF%.bcf}"_tagfilled.bcf