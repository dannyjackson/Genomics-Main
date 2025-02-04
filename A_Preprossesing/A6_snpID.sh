#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script computes average depth statistics of each sample from the output of A4_indelrealignment.sh

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -b  Path to bam directory for analysis (if you are following the full pipeline, this will be ${OUTDIR}/datafiles/indelrealignment/)
  -r  Run name, required for providing a unique name to output files."
    exit 1
fi

# Parse command-line arguments
while getopts p:b:r: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        b) BAMDIR=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

printf "\n\n\n\n"
date
echo "Current script: A2_ClipOverlap.sh"

# Ensure required variables are set
if [ -z "$BAMDIR" ] || [ -z "$ANGSD" ] || [ -z "$SNPPVAL" ] || [ -z "$MINDEPTHIND" ] || [ -z "$MININD" ] || [ -z "$MINQ" ] || [ -z "$MINMAF" ] || [ -z "$MINMAPQ" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi



ls ${BAMDIR}/*bam > ${OUTDIR}/referencelists/${PROJNAME}.bamlist.txt


${ANGSD}/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
    -P 32 -SNP_pval ${SNPPVAL} -setMinDepthInd ${MINDEPTHIND} -minInd ${MININD} -minQ ${MINQ} -minMaf ${MINMAF} -minMapQ ${MINMAPQ} \
    -bam ${OUTDIR}/referencelists/${RUNNAME}.bamlist.txt \
    -out {OUTDIR}/referencelists/${RUNNAME}.allsnps \
    -nThreads ${THREADS}

zcat ${OUTDIR}/referencelists/${RUNNAME}.allsnps.mafs.gz | awk '{print $1, $2, $3, $4}' > ${OUTDIR}/referencelists/${RUNNAME}.sites.mafs

grep ${CHRLEAD} ${OUTDIR}/referencelists/${RUNNAME}.sites.mafs | tail -n +2 > ${OUTDIR}/referencelists/${RUNNAME}.sites_headless.mafs

${ANGSD}/angsd sites index ${OUTDIR}/referencelists/${RUNNAME}.sites_headless.mafs