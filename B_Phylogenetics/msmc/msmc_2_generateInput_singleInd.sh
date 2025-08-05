#!/bin/bash

# Check for at least one argument (parameter file path)
if [ $# -lt 1 ]; then
    echo "Usage: $0 -p <path_to_parameter_file>"
    echo "This script generates a mask for a reference genome to determine which regions are 'snpable'."
    echo "Required Argument:"
    echo "  -p   Path to parameter file (example in GitHub repository as params.sh)"
    echo "  -m   File Name of your unique project msmc params file"
    exit 1
fi

# Parse command-line arguments
while getopts "pm" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        m) MSMCPARAMS=${OPTARG} ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Ensure parameter file is provided and exists
if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
elif [ ! -f "${PARAMS}" ]; then
    echo "Error: Parameter file '${PARAMS}' not found." >&2
    exit 1
fi

# Source the parameter file
source "${PARAMS}"
# Check available modules (useful for debugging environment)
module list

# Source MSMC params file
source "${SCRIPTDIR}/${MSMCPARAMS}"

# Check if filenames list exists
if [ ! -f "${FILENAME_LIST}" ]; then
    echo "Error: Filename list '${FILENAME_LIST}' not found." >&2
    exit 1
fi

# Process each individual sample
while read -r IND; do
    echo "Processing individual: $IND"

    # Check if scaffold list file exists
    if [ ! -f "${OUTDIR}/referencelists/SCAFFOLDS.txt" ]; then
        echo "Error: Scaffold list 'SCAFFOLDS.txt' not found." >&2
        exit 1
    fi

    # Process each scaffold
    while read -r SCAFFOLD; do
        echo "Working on scaffold: $SCAFFOLD"

        VCF="${INDIR}/datafiles/vcf2/${IND}.${SCAFFOLD}.phased.vcf.gz"

        MSMC_INPUT="${MSMCDIR}/input/msmc_input.${IND}.${SCAFFOLD}.txt"

        echo "-----------------------------------------"
        echo "Processing MSMC input generation"
        echo "Individual: ${IND}"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Phasing: ${PHASING}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"
        echo "VCF: ${VCF}"
        echo "-----------------------------------------"

        # Check if the VCF file exists
        if [ -f "$VCF" ]; then
            echo "VCF exists, generating MSMC2 input!"

            # Define mask files
            MASK_INDIV="${MSMCDIR}/mask/ind/ind_mask.${IND}.${SCAFFOLD}.bed.gz"
            MASK_GENOME="${MSMCDIR}/mask/genom/${prefix}_revised_${SCAFFOLD}_mask.${k}.50.bed.gz"

            echo "Individual Mask: ${MASK_INDIV}"
            echo "Genome-wide Mappability Mask: ${MASK_GENOME}"

            if [ "$METHOD" == "samtools" ]; then
                echo "Creating MSMC input file using individual mask (samtools method)."

                ${MSMCTOOLS}/generate_multihetsep.py --mask="${MASK_INDIV}" --mask="${MASK_GENOME}" "${VCF}" > "${MSMC_INPUT}"

            else
                # NOTE: This method was modified on 10 FEB 2024 and may need verification.
                echo "Creating individual mask. Ensure your input VCF includes ALL sites (variant & invariant)."

                VCF_OUT="${VCF}.parsed.vcf"

                # Parse VCF file
                ${MSMCTOOLS}/vcfAllSiteParser.py "$SCAFFOLD" "$MASK_INDIV" "$VCF_OUT"

                echo "Creating MSMC input file with new individual mask."

                ${MSMCTOOLS}/generate_multihetsep.py --mask="${MASK_INDIV}" --mask="${MASK_GENOME}" "${VCF}" > "${MSMC_INPUT}"
            fi

        else
            echo "Warning: VCF file '${VCF}' does not exist, skipping scaffold ${SCAFFOLD}."
        fi

        echo "Finished processing scaffold ${SCAFFOLD}."
    done < ${OUTDIR}/referencelists/SCAFFOLDS.txt

    echo "Finished processing individual ${IND}."

done < "$FILENAME_LIST"

echo "Script execution completed."
date
