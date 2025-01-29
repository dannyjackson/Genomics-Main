#!/bin/sh

# dxy script (1/N)

if [ $# -lt 1 ]
  then
    echo "This script computes dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide dxy and produce the output files necessary for sliding window dxy and dxy for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as fst_params.sh"

  else
    while getopts p:w:s: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;

    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"


printf "\n \n \n \n"
date
echo "Current script: dxy_1.sh"


mkdir -p ${OUTDIR}/analyses/dxy/
mkdir -p ${OUTDIR}/analyses/dxy/${POP1}_${POP2}

cd ${OUTDIR}/analyses/dxy/${POP1}_${POP2}

module load R/4.4.0

total_lines=$(zcat ${OUTDIR}/datasets/safs/${POP1}.mafs.gz | wc -l)
num_sites=$((total_lines - 1))

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p ${OUTDIR}/datasets/safs/${POP1}.mafs.gz -q ${OUTDIR}/datasets/safs/${POP2}.mafs.gz -t ${num_sites}

mv Dxy_persite.txt Dxy_persite_${POP1}_${POP2}.txt

awk '{print $2}' Dxy_persite_${POP1}_${POP2}.txt  | tail -n +2 > ${POP1}_${POP2}_sites.txt

tail -n +2 Dxy_persite_${POP1}_${POP2}.txt > Dxy_persite_${POP1}_${POP2}.forplot.txt

Rscript ${scriptdir}/plotDXY.R -i Dxy_persite_${POP1}_${POP2}.forplot.txt -o ${POP1}_${POP2}_windowed -p ${POP1}_${POP2}_sites.txt -w ${WIN} -s ${STEP}

