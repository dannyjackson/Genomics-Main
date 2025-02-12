#!/bin/sh

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <parameter_file> -i <input_file> -n <population_name> -m <metric> -w <window_size>

This script creates a gene list given a list of sites in a genome and a gff file.

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  Path to the parameter file (e.g., params_genes.sh)
    -i  Path to input file
    -n  Name of population
    -m  Metric (e.g. fst, dxy, Tajima)
    -w  Window size for sliding window analysis 
EOF
    exit 1
fi


# Parse command-line arguments
while getopts "p:i:n:m:w:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) IN_FILE=${OPTARG} ;;
        n) POP=${OPTARG} ;;
        m) METRIC=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        *) echo "Error: Invalid option '-${OPTARG}'" >&2; exit 1 ;;
    esac
done

# Defined in params file:
# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff

# load modules

module load python/3.11/3.11.4
# has to be done on elgato for now
module load gnu8/8.3.0
module load bedtools2/2.29.2

# load parameters
source ${PARAMS}

# defined in generalized
OUT_FILE=${OUTDIR}/analyses/genelist/bed/${POP}.${METRIC}.${WIN}kb.bed
GENES_FILE=${OUTDIR}/analyses/genelist/genes/${POP}.${METRIC}.${WIN}kb.genes.txt
GENENAMES=${OUTDIR}/analyses/genelist/gene_names/${POP}.${METRIC}.${WIN}kb.genenames.txt


python ${SCRIPTDIR}/Genomics-Main/general_scripts/outlier_to_bed.py ${IN_FILE} ${WIN} ${OUT_FILE} ${CHR_FILE}

bedtools intersect -a ${GFF} -b ${OUT_FILE} -wa > ${GENES_FILE}

grep 'gene' ${GENES_FILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}