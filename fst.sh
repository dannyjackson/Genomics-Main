#!/bin/sh

# FST script
# This file accesses fst_window.r and fst_snps.r

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst, sliding window Fst, and fst for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    If you have run it once, feel free to change the window size and rerun. It will skip redundant analyses if the output files already exist (including Fst at individual SNPs) and only compute windowed analyses.

    [-o] Output directory
    [-p] Path to scripts
    [-x] Name of population 1
    [-y] Name of population 2
    [-r] Path to reference genome fasta
    [-g] Path to reference genome gff file

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)
    [-a] Path to angsd directory, will install angsd here if it is not already installed! (defaults to ~/program/angsd/)
    [-c] Leading text before numbers of chromosomes (include only chromosomes not scaffolds; defaults to "NC_0")"

  else
    while getopts i:r:t:p:b:s: option
    do
    case "${option}"
    in
    o) OUTDIR=${OPTARG};;
    p) scriptdir=${OPTARG};;
    x) POP1=${OPTARG};;
    y) POP2=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;
    a) ANGSD=${OPTARG};;
    c) CHRLEAD=${OPTARG};;
    esac
    done

WIN="${WIN:-10000}"
STEP="${STEP:-10000}"
ANGSD="${ANGSD:-~/programs/angsd/}"
CHRLEAD="${CHRLEAD:-NC_0}"

module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2
