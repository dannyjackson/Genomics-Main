#!/bin/bash

# Ensure at least one argument is provided
if [ $# -lt 1 ]; then
    cat <<EOF
Usage: $(basename "$0") -p <population_name> 

This script filters a GFF to just those with enough data across an analysis performed in RAiSD to be considered part of the true background gene list.

Ensure you read and understand the entire script before running it.

REQUIRED ARGUMENTS:
    -p  pop
EOF
    exit 1
fi

# Parse command-line arguments
while getopts "p:" option; do
    case "${option}" in
        p) pop=${OPTARG} ;;
    esac
done

module load bedtools2

GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff # path to gff file

MAPFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.raisd.${pop}.map.filtered.bam"

DEPTHFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.raisd.${pop}.depth.filtered.bam"

MAPBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/raisd.${pop}.map.windows.bed"
DEPTHBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/raisd.${pop}.depth.windows.bed"

BACKGROUNDFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/genes/background_genes.raisd.${pop}.txt"

BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.${pop}.raisd.txt"


tr ' ' '\t' < ${MAPFILE} > ${MAPBED}
tr ' ' '\t' < ${DEPTHFILE} > ${DEPTHBED}

echo "made bedfiles for ${pop}"
echo "intersecting gff and mapbed for ${pop}"

bedtools intersect -a ${GFF} -b ${MAPBED} | sort -u > "${BACKGROUNDFILE}.justmap"

echo "intersecting map-filtered gff and depthbed for ${pop}"

bedtools intersect -a "${BACKGROUNDFILE}.justmap" -b ${DEPTHBED} | sort -u > ${BACKGROUNDFILE}

echo "printing gene list for ${pop}"

grep 'ID\=gene' ${BACKGROUNDFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${BACKGROUNDGENES}