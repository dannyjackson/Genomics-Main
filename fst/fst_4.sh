#!/bin/sh

# FST script (4/4)

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst and produce the output files necessary for sliding window Fst and fst for each SNP. 
    
    Read and understand the entire script before running it!

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as params.sh)

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)"

  else
    while getopts p:w:s: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;

    esac
    done

source ${PARAMS}


printf "\n \n \n \n"
date
echo "Current script: fst_3.sh"



echo "Making relevant gene lists"

# snps 

awk 'BEGIN {FS = ","} {$1=""}1' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv.tmp

mv ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv.tmp ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv

sed -i 's/\"//g' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv

tail -n +2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv  > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.headless.csv


awk -F',' 'NR>1 {print $CHRLEAD, $1, ".1" "\t" $2-1 "\t" $2}' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.headless.csv > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.bed

bedtools intersect -a ${GFF} -b ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.bed -wa > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenes_snps_top.95.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenes_snps_top.95.txt | sed 's/ID\=gene\-//g' | sort -u > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenenames_snps_top.95.txt

# windowed

awk 'BEGIN {FS = ","} {$1=""}1' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv.tmp

mv ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv.tmp ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv

sed -i 's/\"//g' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv

tail -n +2 ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.headless.csv


awk -F',' -v win="$WIN" 'NR>1 {print $CHRLEAD, $1, ".1" "\t" $2-(win/2) "\t" $2+(win/2)}' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.headless.csv > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.outlierfst.bed

bedtools intersect -a ${GFF} -b ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.outlierfst.bed -wa > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenes_windowed_top.95.txt

# awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenes_windowed_top.95.txt | sed 's/ID\=gene\-//g' | sort -u > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenenames_windowed_top.95.txt

fi