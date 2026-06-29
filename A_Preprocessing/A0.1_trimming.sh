#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script trims FASTQ files with options to use various trimming methods.

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).

Notes on trimming methods:
    - trimmomatic: Ensure you have the trimmomatic jar in programs directory
    - trimgalore & fastp: This script assumes that you have installed and loaded these programs using a conda/micromamba environmenty 

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -i Individual name (can easily be passed through a slurm array).
  -m Optional Argument to state method of trimming. Options include "trimmomatic", "trimgalore", "fastp" (Defaults to trimmomatic)"
    exit 1
fi

METHOD=trimmomatic

# Parse command-line arguments
while getopts p:i:m: option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
		i) IND=${OPTARG};;
        m) METHOD=${OPTARG};;
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
echo "Current script: A0.1_trimming.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$FASTAS" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


if [ "$METHOD" = "trimmomatic" ]; then

    echo "Beginning trimming for "$IND "using trimmomatic">>${OUTDIR}/datafiles/trimming/${IND}_trim_log.txt
    java -jar ${TRIMJAR} PE -threads ${THREADS} ${FASTAS}/"$IND"_R1.fastq.gz  ${FASTAS}/"$IND"_R2.fastq.gz  \
    -baseout ${OUTDIR}/datafiles/trimmed_fastas/"$IND"_trimmed.fq.gz \
    LEADING:${LEAD} TRAILING:${TRAIL} SLIDINGWINDOW:${SLIDE} MINLEN:${MINREADLEN} >> ${OUTDIR}/datafiles/trimming/${IND}_trim_log.txt

elif [ "$METHOD" = "trimgalore" ]; then
    echo "Beginning trimming for "$IND "using trimgalore">>${OUTDIR}/datafiles/trimming/${IND}_trim_log.txt
    trim_galore --paired --cores ${THREADS} ${FASTAS}/${IND}_R1.fastq.gz ${FASTAS}/${IND}_R2.fastq.gz -o ${OUTDIR}/datafiles/trimmed_fastas
    # Renaming output files to be consistent with pipeline
    mv ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R1_val_1.fq.gz ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R1_trimmed.fq.gz
    mv ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R2_val_2.fq.gz ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R2_trimmed.fq.gz

elif [ "$METHOD" = "fastp" ]; then

    echo "Beginning trimming for "$IND "using fastp">>${OUTDIR}/datafiles/trimming/${IND}_trim_log.txt
    fastp --in1 ${FASTAS}/"${IND}"_R1.fastq.gz --in2 ${FASTAS}/"${IND}"_R2.fastq.gz -q 5 -n 15 \
    --out1 ${OUTDIR}/datafiles/trimmed_fastas/"${IND}"_R1_trimmed.fq.gz --out2 ${OUTDIR}/datafiles/trimmed_fastas/"${IND}"_R2_trimmed.fq.gz

fi

echo Finished trimming ${IND} using ${METHOD}
