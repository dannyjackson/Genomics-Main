#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script computes average depth statistics of a bam file. Meant to be used in a slurm array

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -r  Run name, required for providing a unique name to output files.
  -i  Individual name (can easily be passed through a slurm array).
  -c  Optional flag to indicate whether to compress the output files with bgzip."
    exit 1
fi

# Set default values for optional variables
COMPRESS=false

# Parse command-line arguments
while getopts p:r:i:c option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        i) IND=${OPTARG};;
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
echo "Current script: A0.5_bamstatistics_array.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] ; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi

## Compute statistics on bam file
echo "Current bam: $IND"
echo $IND >> ${OUTDIR}/datafiles/bamstats/"$IND"_depthstats.txt 
samtools depth ${BAMDIR}/"$IND".realigned.bam >> ${OUTDIR}/datafiles/bamstats/"$IND"_depthstats.txt 
echo "Finished: $IND"


# Save summary stats in CSV format
echo "Generating summary bamstats for $IND"
# Compute average and standard deviation
stats=$(awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR "," sqrt(sumsq/NR - (sum/NR)**2) }' ${OUTDIR}/datafiles/bamstats/"$IND"_depthstats.txt)

if [[ -f "${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt" ]]; then
    echo "$IND,$stats" >> ${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt
else
    echo "Sample,Average,Stdev" > ${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt
    echo "$IND,$stats" >> ${OUTDIR}/datafiles/bamstats/${RUNNAME}.depthstats.txt
fi


# Since the sample bamstats files are really large, let's compress them for now while we may need them in the near future.
if [ "$COMPRESS" = true ]; then
    echo "Compressing ${IND} depthstats file"
    bgzip "${OUTDIR}/datafiles/bamstats/${IND}_depthstats.txt"
fi