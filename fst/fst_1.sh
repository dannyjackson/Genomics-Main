#!/bin/sh

# FST script (1/4)

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the A8_siteallelefrequency.sh script in github.com/dannyjackson/Genomics-Main/A_Preprocessing. It will compute average genome-wide Fst and produce the output files necessary for sliding window Fst and Fst for each SNP. 
    
    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in the GitHub repository as params.sh)."
  else
    while getopts p: option
    do
    case "${option}" in
    p) PARAMS=${OPTARG};;
    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"

printf "\n \n \n \n"
date
echo "Current script: fst.sh"

if [ -f "${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml" ]
        then
            echo "Files present in MLS directory, assuming they are already generated and moving on!"
        else
            echo "Computing the 2D SFS prior"
            ${ANGSD}/misc/realSFS ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx > ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml
fi

if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}"* ]
        then
            echo "FST index file already present, assuming it is already generated and moving on!"
        else
            echo "Preparing FST for window analysis, computing FST index file"
            ${ANGSD}/misc/realSFS fst index ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx -sfs ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml -fstout ${OUTDIR}/analyses/fst/${POP1}_${POP2}
fi

if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt" ]
        then
            echo "Global FST estimate file already present, assuming it is already generated and moving on!"
        else
            echo "Computing global FST estimate"
            echo -e "FST.Unweight\tFST.Weight" > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt
            ${ANGSD}/misc/realSFS fst stats ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx >> ${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt
fi

fi
