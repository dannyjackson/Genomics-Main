#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script creates a masking file and a vcf file for a given population."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -b  List of BAMs to process."
    echo "  -v  Path to output VCF file."
    exit 1
}

# Parse command-line arguments
while getopts ":p:b:v:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        b) BAMLIST=${OPTARG} ;;
        v) VCF_OUT=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Validate required arguments
if [[ -z "${PARAMS:-}" || -z "${BAMLIST:-}" || -z "${VCF_OUT:-}" ]]; then
  echo "Error: Missing required arguments." >&2
  usage
fi

if [[ ! -f "${PARAMS}" ]]; then
  echo "Error: Parameter file '${PARAMS}' not found." >&2
  exit 1
fi
if [[ ! -f "${BAMLIST}" ]]; then
  echo "Error: BAM list file '${BAMLIST}' not found." >&2
  exit 1
fi

# Load parameters
source "${PARAMS}"


# THREADS and REF must be defined in params_preprocessing.sh
: "${REF:?REF not set in params file}"
: "${THREADS:?THREADS not set in params file}"

# Print script info
echo
date
echo "Params:   ${PARAMS}"
echo "REF:      ${REF}"
echo "THREADS:  ${THREADS}"
echo "BAMLIST:  ${BAMLIST}"
echo "VCF_OUT:  ${VCF_OUT}"
echo

# Check BAMs listed exist + index if needed
while read -r BAM; do
  [[ -z "${BAM}" ]] && continue
  [[ "${BAM}" =~ ^# ]] && continue
  if [[ ! -f "${BAM}" ]]; then
    echo "Error: BAM not found: ${BAM}" >&2
    exit 1
  fi
  if [[ ! -f "${BAM}.bai" ]]; then
    echo "Indexing: ${BAM}"
    samtools index "${BAM}"
  fi
done < "${BAMLIST}"

# Ensure reference indexed
if [[ ! -f "${REF}.fai" ]]; then
  echo "Indexing reference fasta..."
  samtools faidx "${REF}"
fi

# Generate VCF: biallelic SNPs only, bgzipped output
bcftools mpileup -Ou -f "${REF}" -b "${BAMLIST}" --threads "${THREADS}" \
  | bcftools call -mv -Ou --threads "${THREADS}" \
  | bcftools view -v snps -m2 -M2 -Oz -o "${VCF_OUT}"

tabix -p vcf "${VCF_OUT}"

# Final report
echo "Done. Wrote: ${VCF_OUT}"
date