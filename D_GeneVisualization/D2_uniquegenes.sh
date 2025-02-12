#!/bin/sh

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <parameter_file> -i <pop1_name> -q <pop2_name> -m <metric> -w <window_size>

This script creates a gene list given a list of sites in a genome and a GFF file.

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  Path to parameter file
    -i  Name of the population of interest
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
        *) echo "Error: Invalid option '-${option}'" >&2; exit 1 ;;
    esac
done

# Ensure the parameter file exists
if [ ! -f "$PARAMS" ]; then
    echo "Error: Parameter file '$PARAMS' not found!" >&2
    exit 1
fi

source "$PARAMS"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not set. Ensure the parameter file defines OUTDIR." >&2
    exit 1
fi

# Define file paths
GENE_LIST1="$OUTDIR/analyses/genelist/gene_names/${POP1}.${METRIC}.${WIN}.genenames.txt"
GENE_LIST2="$OUTDIR/analyses/genelist/gene_names/${POP2}.${METRIC}.${WIN}.genenames.txt"
FINAL_LIST="$OUTDIR/analyses/genelist/final_gene_lists/${POP1}.${METRIC}.${WIN}.unique.genenames.txt"

# Ensure gene list files exist
if [ ! -f "$GENE_LIST1" ]; then
    echo "Error: Gene list for $POP1 not found: $GENE_LIST1" >&2
    exit 1
fi

if [ ! -f "$GENE_LIST2" ]; then
    echo "Warning: Gene list for $POP2 not found. Keeping all genes from $POP1."
    cp "$GENE_LIST1" "$FINAL_LIST"
else
    grep -Fxv -f "$GENE_LIST2" "$GENE_LIST1" > "$FINAL_LIST"
fi
