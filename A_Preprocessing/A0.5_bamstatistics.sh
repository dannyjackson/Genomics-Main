#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script computes average depth statistics of each bam file in a directory.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -r  Run name, required for providing a unique name to output files.
  -c  Optional flag to indicate whether to compress the output files with bgzip."
    exit 1
fi

# Set default values for optional variables
COMPRESS=false

# Parse command-line arguments
while getopts p:r:c option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        c) COMPRESS=true ;;
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
echo "Current script: A0.5_bamstatistics.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] ; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


## Create a list of sample ids for this run assumes bams are all formated with sample names as first term in an underscore separated string
ls ${BAMDIR} | awk -F "." '{print $1}' | sort -u > "${OUTDIR}/referencelists/${RUNNAME}.sampleids.txt"


## Compute statistics on bam files 

while read -r bird; do
echo "Current bam: $bird"
echo $bird >> ${OUTDIR}/datafiles/bamstats/"$bird"_depthstats.txt 
samtools depth ${BAMDIR}/"$bird".realigned.bam >> ${OUTDIR}/datafiles/bamstats/"$bird"_depthstats.txt 
echo "Finished: $bird"
done <  ${OUTDIR}/referencelists/${RUNNAME}.sampleids.txt


# Create one file with all bam stats
echo "Generating summary bamstats file for all individuals"

# Print header
echo "Sample,Average,Stdev" > ${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt

while read -r bird; do 
  # Compute average and standard deviation
  stats=$(awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR "," sqrt(sumsq/NR - (sum/NR)**2) }' ${OUTDIR}/datafiles/bamstats/"$bird"_depthstats.txt)
  
  # Append results in CSV format
  echo "$bird,$stats" >> ${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt

done < ${OUTDIR}/referencelists/${RUNNAME}.sampleids.txt


# Since the sample bamstats files are really large, let's compress them for now while we may need them in the near future.
if [ "$COMPRESS" = true ]; then
    for bird in $(cat ${OUTDIR}/referencelists/${RUNNAME}.sampleids.txt); do
         echo "Compressing ${bird} depthstats file"
         bgzip "${OUTDIR}/datafiles/bamstats/${bird}_depthstats.txt"
    done
fi