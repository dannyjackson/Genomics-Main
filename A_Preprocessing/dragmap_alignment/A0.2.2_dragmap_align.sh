#!/bin/sh

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script aligns FASTQ files to a reference genome using Dragmap Open Source Aligner (https://github.com/Illumina/dragmap).

I recommend running it as a slurm array to pass individuals to sbatch jobs for maximum efficiency (see github.com/dannyjackson/BioinformaticTutorials/SubmittingJobs.txt for an explanation of running slurm arrays).
This script assumes that you've installed dragmap using conda/micromamba environment

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -i  Individual name (can easily be passed through a slurm array)."
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
echo "Current script: A0.2.2_dragmap_align.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$FASTAS" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi

echo "Making hashtable from ${REF}"
mkdir ${OUTDIR}/referencelists/dragmap_reference
dragen-os --build-hash-table true --ht-reference ${REF}  --output-directory ${OUTDIR}/referencelists/dragmap_reference
echo "Hash table stored in ${OUTDIR}/referencelists/dragmap_reference"

echo "Aligning ${IND} with dragmap"
dragen-os -r ${OUTDIR}/referencelists/dragmap_reference/ -1 ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R1_trimmed.fq.gz -2 ${OUTDIR}/datafiles/trimmed_fastas/${IND}_R2_trimmed.fq.gz \
    --output-directory ${OUTDIR}/datafiles/bamfiles/  --output-file-prefix ${IND}

echo "Done aligning ${IND}"

echo "Converting Dragmap SAM to BAM"
samtools view -bS ${OUTDIR}/datafiles/bamfiles/${IND}.sam > ${OUTDIR}/datafiles/bamfiles/${IND}.bam
rm ${OUTDIR}/datafiles/bamfiles/${IND}.sam
echo "Deleted Intermediate File ${IND}.sam"

echo "Finished ${IND}"