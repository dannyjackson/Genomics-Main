#!/bin/sh

# FST script
# This file accesses fst_window.r and fst_snps.r

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst and produce the output files necessary for sliding window Fst and fst for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as fst_params.sh)

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)

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

source ${PARAMS}

WIN="${WIN:-10000}"
STEP="${STEP:-10000}"


printf "\n \n \n \n"
date
echo "Current script: fst.sh"

# Generate saf files for each population in angsd. Skip if there is any output already in ${OUTDIR}/datafiles/safs/.

if [ -f "${OUTDIR}/datafiles/safs/" ]
        then
            echo "Files present in safs directory, assuming they are already generated and moving on!"
        else
            echo "computing safs for pop 1"
            ~/programs/angsd/angsd -bam ${OUTDIR}/datafiles/${POP1}.bams.txt -out ${OUTDIR}/datafiles/safs/${POP1} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites ${OUTDIR}/referencelists/sites_headless.mafs -anc ${REF} -nThreads 10
            
            echo "computing safs for pop 2"
            ~/programs/angsd/angsd -bam ${OUTDIR}/datafiles/${POP2}.bams.txt -out ${OUTDIR}/datafiles/safs/${POP1} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites ${OUTDIR}/referencelists/sites_headless.mafs -anc ${REF} -nThreads 10
fi


if [ -f "${OUTDIR}/datafiles/mls/" ]
        then
            echo "Files present in mls directory, assuming they are already generated and moving on!"
        else
            echo "computing the 2dsfs prior"
            ${ANGSD}/misc/realSFS ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx > ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml

fi


if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}*" ]
        then
            echo "fst.idx file already present, assuming it is already generated and moving on!"
        else
            echo "prepare the fst for easy window analysis, compute fst.idx file"
            ${ANGSD}/misc/realSFS fst index ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx -sfs ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml -fstout ${OUTDIR}/analyses/fst/${POP1}_${POP2}
fi

if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt" ]
        then
            echo "global fst estimate file already present, assuming it is already generated and moving on!"
        else
            echo "computing global fst estimate"
            ${ANGSD}/misc/realSFS fst stats ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt

fi
