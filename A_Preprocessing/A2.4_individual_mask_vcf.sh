#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script computes average depth statistics of each BAM file in a directory."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the individual to analyze."
    echo "  -b  Directory containing BAM files."
    exit 1
}

# Parse command-line arguments
while getopts ":p:i:b:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) IND=${OPTARG} ;;
        b) BAMDIR=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Validate required arguments
if [[ -z "${PARAMS}" || -z "${IND}" || -z "${BAMDIR}" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

# Ensure the parameter file exists
if [[ ! -f "${PARAMS}" ]]; then
    echo "Error: Parameter file '${PARAMS}' not found." >&2
    exit 1
fi

# Load parameters
source "${PARAMS}"

# Define BAM file
BAMFILE="${BAMDIR}/${IND}.realigned.bam"

# Print script info
echo -e "\n\n"
date
echo "Individual: ${IND}"
echo "BAM file: ${BAMFILE}"

# Check and create BAM index if needed
if [[ -f "${BAMFILE}.bai" ]]; then
    echo "BAM index found, proceeding..."
else
    echo "BAM index not found, creating index..."
    samtools index "${BAMFILE}"
fi

# Ensure scaffold list file exists
SCAFFOLD_FILE="${OUTDIR}/referencelists/SCAFFOLDS.txt"
if [[ ! -f "${SCAFFOLD_FILE}" ]]; then
    echo "Error: Scaffold list file '${SCAFFOLD_FILE}' not found." >&2
    exit 1
fi

echo "processing scaffolds" 
# Process each scaffold
while read -r s; do
    echo "Processing scaffold: ${s}"

    # Calculate mean coverage
    MEANCOV=$(samtools depth -r "${s}" "${BAMFILE}" | 
              awk '{sum += $3} END {if (NR==0) print 0; else print sum / NR}' | tr ',' '.')

    echo "${IND}.${s} ${MEANCOV}" >> "${OUTDIR}/datafiles/bamstats/coverage_samtoolsDepth_${IND}.txt"
    echo "Mean coverage for scaffold ${s}: ${MEANCOV}"

    # Define output file paths
    MASK_IND="${OUTDIR}/datafiles/mask/ind_mask.${IND}.${s}.bed.gz"
    VCF="${OUTDIR}/datafiles/vcf/${IND}.${s}.vcf"

    # Ensure reference genome is indexed
    if [[ ! -f "${REF}.fai" ]]; then
        echo "Indexing reference genome..."
        samtools faidx "${REF}"
    fi

    # Generate VCF and mask file if using samtools

    echo "Running samtools for scaffold ${s}..."
    bcftools mpileup -Ou -r "${s}" --threads "${THREADS}" -f "${REF}" "${BAMFILE}" | \
    bcftools call -mv -V indels --threads "${THREADS}" | \
    "${SCRIPTDIR}/msmc_tools/bamCaller.py" "${MEANCOV}" "${MASK_IND}" > "${VCF}"

    echo "Completed scaffold ${s}."


done < "${SCAFFOLD_FILE}" 


# Final report
echo "Filtered VCF and mask created for all scaffolds for ${IND}."
echo "Script completed."
date
