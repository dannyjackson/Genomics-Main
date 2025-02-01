#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file>

This script realigns around indels using GATK 3.7.0.

GATK retired the indelrealignment command in later versions, and to my knowledge no other resource exists for this function. GATK just incorporated the indel realignment into their align-call pipeline, so if you use GATK for SNP calling, it is just fine. However, for angsd pipelines, we have to rely on this old version of GATK for indel realignment.

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
echo "Current script: A2_ClipOverlap.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$FASTAS" ] || [ -z "$BAMUTILBAM" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi


# Create necessary directories
mkdir -p "${OUTDIR}/datafiles/indelmaps"



# Index bams
samtools index ${OUTDIR}/datafiles/clipoverlap/$IND.all.sorted.marked.clipped.bam 
echo "done " ${IND} >>${OUTDIR}/datafiles/clipoverlap/index_clippedstats.txt 

# Create indel maps
apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R ${REF} \
-I ${OUTDIR}/datafiles/clipoverlap/${IND}.all.sorted.marked.clipped.bam \
-o ${OUTDIR}/datafiles/indelmaps/${IND}.intervals

# Realign around indels
apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
-R ${REF} \
--consensusDeterminationModel USE_READS \
-I ${OUTDIR}/datafiles/clipoverlap/${IND}.all.sorted.marked.clipped.bam \
--targetIntervals ${OUTDIR}/datafiles/indelmaps/${IND}.intervals \
-o ${OUTDIR}/datafiles/indelrealignment/${IND}.realigned.bam
