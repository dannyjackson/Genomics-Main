#!/bin/bash

# Ensure script directory and output directory are defined
if [ -z "$SCRIPTDIR" ]; then
    echo "Error: SCRIPTDIR is not defined. Please set this variable."
    exit 1
fi

# Load parameters
if [ ! -f "${SCRIPTDIR}/params.sh" ]; then
    echo "Error: params.sh not found in ${SCRIPTDIR}"
    exit 1
fi
source "${SCRIPTDIR}/params.sh"


# Process each individual sample
while read -r IND; do
    echo "Processing individual: $IND"
    BAMFILE="${BAMDIR}/${IND}.final.bam"
    if [[ ! -f "$BAMFILE" ]]; then
        echo "Warning: BAM file for $IND not found. Skipping."
        continue
    fi

    # Iterate through scaffolds
    while read -r SCAFFOLD; do
        echo "Processing scaffold: $SCAFFOLD"
        VCF_OUT="${OUTDIR}/vcf/${IND}.${SCAFFOLD}.${phasing}.samtools.vcf.gz"
        VCF_IN="${OUTDIR}/vcf2/${IND}.${SCAFFOLD}.samtools.vcf"

        if [[ -f "$VCF_OUT" ]]; then
            echo "Phased VCF already exists for scaffold $SCAFFOLD; skipping."
        else
            echo "Phasing VCF for scaffold $SCAFFOLD..."
            sed -i.bak 's/^ //g' "$VCF_IN"
            
            whatshap phase --reference "$GENOME" --ignore-read-groups \
                          -o "$VCF_OUT" "$VCF_IN" "$BAMFILE"
            
            whatshap stats --tsv="${OUTDIR}/stats/${IND}.${SCAFFOLD}.${prefix}.minDP10.${phasing}.stats.tsv" \
                          "$VCF_OUT"
        fi
    done < "${OUTDIR}/referencelists/SCAFFOLDS.txt"

done < /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt

echo "VCF phasing completed."
date
