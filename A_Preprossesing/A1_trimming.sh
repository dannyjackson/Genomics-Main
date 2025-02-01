#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script aligns FASTQ files to a reference genome using BWA MEM.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -i Individual name (can easily be passed through a slurm array)."
    exit 1
fi

# Parse command-line arguments
while getopts p:i: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		i) IND=${OPTARG};;
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
echo "Current script: A2_alignandsort.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$FASTAS" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/condensed_fastas" "${OUTDIR}/datafiles/trimming"
bwa index ${REF}

# create list of samples, assumes fastas are all formated with sample names as first term in an underscore separated string
ls ${FASTAS} | awk -F "_" '{print $1}' | sort -u > "${OUTDIR}/referencelists/sampleids.txt"


while read -r ID;
do
  echo "Beginning trimming for "$ID>>${OUTDIR}/trimming/trim_log.txt
  java -jar ${TRIMJAR} PE -threads 12 \
  ${OUTDIR}/condensed_fastas/"$ID"_1.fq.gz  ${OUTDIR}/condensed_fastas/"$ID"_2.fq.gz  \
  -baseout ${OUTDIR}/datafiles/trimmed_fastas/"$ID"_trimmed.fq.gz \
  LEADING:${LEAD} TRAILING:${TRAIL} SLIDINGWINDOW:${SLIDE} MINLEN:${MINREADLEN}>>${OUTDIR}/datafiles/trimming/trim_log.txt


done < "${OUTDIR}/referencelists/sampleids.txt"
