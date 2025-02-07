#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> -b <bam_directory>"
    echo ""
    echo "This script combines all vcfs from a single scaffold, where each individual has its own vcf, into a single fasta file. For use in phylonetwork and other phylogenetic analyses"
    echo "It is best run as a Slurm array that calls this script for each scaffold."
    echo ""
    echo "Required arguments:"
    echo "  -s  Name of the scaffold to analyze."
    exit 1
}


# Parse command-line arguments
while getopts ":s:" option; do
    case "${option}" in
        s) SCAFFOLD=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done


REF="/xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna"
IND_LIST="/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt"

# Validate that the individual list file exists
if [[ ! -f "$IND_LIST" ]]; then
    echo "Error: Individual list file not found: $IND_LIST"
    exit 1
fi

echo "Processing scaffold: $SCAFFOLD"
# process each individual
while read -r IND; do
    echo "Processing individual: $IND"
    
    # Define input VCF file path
    VCF_IN="${OUTDIR}/datafiles/vcf2/${IND}.${SCAFFOLD}.phased.vcf.gz"
    
    # Check if the VCF file exists
    if [[ ! -f "$VCF_IN" ]]; then
        echo "Error: VCF file not found: $VCF_IN"
        continue
    fi

    samtools faidx $REF $SCAFFOLD| bcftools consensus /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/${IND}.${SCAFFOLD}..samtools.vcf.gz -s $IND | sed 1d >> /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/"$SCAFFOLD".fa
done < "$IND_LIST"
