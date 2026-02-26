#!/bin/sh

# FST script (1/4)

if [ $# -lt 1 ]
  then
    echo "This script computes the genome-wide site allele frequency spectrum within a population of genomes using a genotype likelihood framework implemented in ANGSD. It requires a list of pre-processed bam files, the scripts for which are A1-A4 in this repository. 
    You should make your population bam list on your own and save them in the format POP_bams.txt in the referencelists directory.

    This is an essential step for analyzing fst, dxy, and Tajima's D.
    
    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in the GitHub repository as params.sh).
    -n Name of population 
    -m Path to headless MAF file generated in A1.1_snpID.sh"

  else
    while getopts p:n:m: option
    do
    case "${option}" in
    p) PARAMS=${OPTARG};;
    n) POP=${OPTARG};;
    m) MAFFILE=${OPTARG};;
    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"

printf "\n \n \n \n"
date
echo "Current script: A1.3_siteallelefrequency.sh"


# Generate SAF files for each population in ANGSD. Skip if output already exists.

if [ -f "${OUTDIR}/datafiles/safs/${POP}"* ]
        then
            echo "${POP} files present in SAFs directory, assuming they are already generated and moving on!"
        else
            echo "Computing SAFs for ${POP}"
            ${ANGSD}/angsd -bam ${OUTDIR}/referencelists/${POP}_bams.txt -out ${OUTDIR}/datafiles/safs/${POP} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 \
            -setMinDepthInd ${MINDEPTHIND} -minInd ${MININD} -minQ ${MINQ} -minMapQ ${MINMAPQ} -sites ${MAFFILE} -anc ${REF} -nThreads ${THREADS} 
fi


fi