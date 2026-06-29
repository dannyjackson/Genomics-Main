#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script computes average depth statistics of a bam file. Meant to be used in a slurm array

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -r  Run name, required for providing a unique name to output files.
  -i  Individual name (can easily be passed through a slurm array).
  -d  Optional flag to indicate whether to delete depthstats file (recommended)."
    exit 1
fi

# Set default values for optional variables
DELETE=false

# Parse command-line arguments
while getopts p:r:i:d option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        i) IND=${OPTARG};;
        d) DELETE=true ;;
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
echo "Current bampath: ${BAMDIR}/$IND.realigned.bam"
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


if [ "$DELETE" = true ]; then
    rm "${OUTDIR}/datafiles/bamstats/${IND}_depthstats.txt"
    echo "Deleted ${IND} depthstats file"
fi

echo "Finished bam statistics for $IND"