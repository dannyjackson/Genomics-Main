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
mkdir -p "${OUTDIR}/datafiles/trimmed_fastas" "${OUTDIR}/datafiles/bamfiles/${IND}" "${OUTDIR}/datafiles/sortedbamfiles/${IND}"

# Align reads using BWA MEM
bwa mem -t "${THREADS}" "${REF}" \
    "${FASTAS}/${IND}_trimmed_1P.fq.gz" \
    "${FASTAS}/${IND}_trimmed_2P.fq.gz" | \
    samtools view -b -o "${OUTDIR}/datafiles/bamfiles/${IND}.bam" -S

# Sort BAM file
samtools sort -T "${OUTDIR}/datafiles/sortedbamfiles/${IND}/temp" -@ "${THREADS}" \
    -o "${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted.bam" \
    "${OUTDIR}/datafiles/bamfiles/${IND}.bam"

# Add read groups
picard AddOrReplaceReadGroups \
    I="${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted.bam" \
    O="${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded.bam" \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${IND}"

# Mark duplicates
picard MarkDuplicates \
    I="${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded.bam" \
    O="${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam" \
    M="${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}.duplicate.metrics.txt"

# Index the final BAM file
samtools index "${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam"

# Generate alignment statistics
samtools flagstat "${OUTDIR}/datafiles/sortedbamfiles/${IND}/${IND}_sorted_RGadded.bam"
