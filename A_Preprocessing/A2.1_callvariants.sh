#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <parameter_file> -b <bamlist> -r <run_name>

This script calls variants for VCF file creation from BAM files.

Required argument:
  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository).
  -b  Path to bam list file for analysis
  -r  Run name, required for providing a unique name to output files.
Optional argument:
  -s  Boolean flag to enable split-region variant calling. Will require two extra files (named chromosomes-id file, and unplaced-scaffolds file generated in base_setup.sh)"
    exit 1
fi

SPLITREGIONS=false

# Parse command-line arguments
while getopts p:r:b:s option; do
    case "${option}" in
        p) PARAMS=${OPTARG};;
        r) RUNNAME=${OPTARG};;
        b) BAMLIST=${OPTARG};;
        s) SPLITREGIONS=true ;;
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
echo "Current script: A2.1_callvariants.sh"

# Ensure required variables are set
if [ -z "$OUTDIR" ] || [ -z "$REF" ] || [ -z "$BAMLIST" ]; then
    echo "Error: Missing required parameters in the parameter file." >&2
    exit 1
fi

if [ "$SPLITREGIONS" = false ]
    then
        echo "Calling variants in single job"
        bcftools mpileup -Ou -f ${REF} -a FORMAT/AD,DP,INFO/AD,SP --bam-list ${BAMLIST} | bcftools call -mv -V indels -o ${OUTDIR}/datafiles/genotype_calls/"${RUNNAME}"_snps_multiallelic.vcf


    else
        echo "Calling variants in parallel across chromosome regions"

        module load parallel
        CHROMS=${OUTDIR}/referencelists/SCAFFOLDS.chroms.txt
        SCAFFOLDS=${OUTDIR}/referencelists/SCAFFOLDS.unlabeled.txt

        # Check for required Chrom/Scaffold files
        if [ -f "${CHROMS}" ]
        then
            echo "List of chromosomes found, continuing..."
        else
            echo "List of chromosomes not found at path: ${CHROMS}. Please fill this file with desired named macrochromosomes (base_setup.sh). EXITING..."
            exit 1
        fi

        if [ -f "${SCAFFOLDS}" ]
        then
            echo "List of unplaced scaffolds found, continuing..."
        else
            echo "List of unplaced scaffolds not found at path: ${SCAFFOLDS}. Please fill this file with desired scaffolds (base_setup.sh). EXITING..."
            exit 1
        fi

        echo "Parallel variant calling on named macrochromosomes..."
        echo "Processing ${THREADS} chromosomes at a time..."
        cat ${CHROMS} | parallel --jobs "$THREADS" "
            bcftools mpileup -f ${REF} -r {} -b ${BAMLIST} -Ou -a FORMAT/AD,DP,INFO/AD,SP | \
            bcftools call -mv -Ou -V indels -o ${OUTDIR}/datafiles/genotype_calls/${RUNNAME}_snps_multiallelic_{}.vcf"

        echo "Calling variants for unplaced scaffolds..."
        bcftools mpileup -f ${REF} -R ${SCAFFOLDS} -b ${BAMLIST} -Ou -a FORMAT/AD,DP,INFO/AD,SP | \
        bcftools call -mv -Ou -V indels -o ${OUTDIR}/datafiles/genotype_calls/${RUNNAME}_snps_multiallelic_scaffolds.vcf

        # concat vcf files in required .fai order
        #bcftools concat -f vcf_list_order.txt -Ou | \
        #bcftools norm -f "$REF" -Ou -o "${OUTDIR}/datafiles/genotype_calls/${RUNNAME}_snps_multiallelic_merged.vcf"

fi

echo "Done"