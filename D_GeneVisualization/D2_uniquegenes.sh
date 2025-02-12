#!/bin/sh

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <pop1_name> -q <pop2_name> -m <metric> -w <window_size>

This script creates a gene list given a list of sites in a genome and a gff file.

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  Name of the population of interest
    -q  Name of the control population (removes genes identified in this pop from the final file)
    -m  Metric (e.g. Tajima, raisd)
    -w  Window size 
EOF
    exit 1
fi


# Parse command-line arguments
while getopts "p:i:q:m:w:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) POP1=${OPTARG} ;;
        q) POP2=${OPTARG} ;;
        m) METRIC=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        *) echo "Error: Invalid option '-${OPTARG}'" >&2; exit 1 ;;
    esac
done

source ${PARAMS}
grep -Fxv -f ${OUTDIR}/analyses/genelist/gene_names/${POP2}.${METRIC}.${WIN}.genenames.txt ${OUTDIR}/analyses/genelist/gene_names/${POP1}.${METRIC}.${WIN}.genenames.txt > ${OUTDIR}/analyses/genelist/final_gene_lists/${POP1}.${METRIC}.${WIN}.unique.genenames.txt
