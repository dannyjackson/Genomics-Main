#!/bin/bash

module load samtools
module load bcftools

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script combines all vcfs from a single scaffold, where each individual has its own vcf, into a single fasta file."
    echo "It is best run as a Slurm array that calls this script for each scaffold."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to raxml params file."
    exit 1
}


# Parse command-line arguments
while getopts ":s:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done


# Validate that the individual list file exists
if [[ ! -f "$IND_LIST" ]]; then
    echo "Error: Individual list file not found: $IND_LIST"
    exit 1
fi

    while read -r SCAFFOLD; do
    echo "Processing scaffold: $SCAFFOLD"
    # process each individual
    while read -r IND; do
        echo "Processing individual: $IND"
        
        # Define input VCF file path
        VCF_IN="/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/${IND}.${SCAFFOLD}..samtools.vcf.gz"
    
        
        # Check if the VCF file exists
        if [[ ! -f "$VCF_IN" ]]; then
            echo "Error: VCF file not found: $VCF_IN"
            continue
        fi
        echo "indexing: {IND}.${SCAFFOLD}..samtools.vcf.gz"
        bcftools index /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/${IND}.${SCAFFOLD}..samtools.vcf.gz

        echo '>' $IND >> /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/"$SCAFFOLD".fa

        echo "Making fasta for: ${IND}.${SCAFFOLD}..samtools.vcf.gz"
        samtools faidx $REF $SCAFFOLD| bcftools consensus /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/${IND}.${SCAFFOLD}..samtools.vcf.gz | sed 1d >> /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/"$SCAFFOLD".fa

    done < "$IND_LIST"
done < "$SCAFFOLD_LIST"
