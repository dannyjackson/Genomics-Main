#!/bin/bash

# Usage message function
usage() {
    echo "Usage: $0 -p <parameter_file> -i <individual> [-s]"
    echo ""
    echo "This script phases a VCF file with SHAPEIT5 using a VCF split by scaffold generated previously through the SNPable pipeline. It can also polish the phased output with SAPPHIRE, which is highly recommended"
    echo "It is best run as a Slurm array that calls this script for each individual."
    echo ""
    echo "Required arguments:"
    echo "  -p  Path to the parameter file (e.g., params_preprocessing.sh in the GitHub repository)."
    echo "  -i  Name of the population to analyze."
    echo "  -s  Optional flag to enable SAPPHIRE polishing of phased VCFs (recommended)."
    exit 1
}

SAPPHIRE=false

# Parse command-line arguments
while getopts ":p:i:s" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        i) POP=${OPTARG} ;;
        s) SAPPHIRE=true ;; # False by default, recommended to set to true
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z "$PARAMS" || -z "$IND" || -z "$BAMDIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Load parameters
source "${PARAMS}"

# Ensure OUTDIR is set
if [ -z "$OUTDIR" ]; then
    echo "Error: OUTDIR is not defined. Please set this variable."
    exit 1
fi

echo "Processing individual: $IND"


# Ensure the scaffold list exists
SCAFFOLD_LIST="${OUTDIR}/referencelists/SCAFFOLDS.txt"
if [[ ! -f "$SCAFFOLD_LIST" ]]; then
    echo "Error: Scaffold list file not found: $SCAFFOLD_LIST"
    exit 1
fi

# Iterate through scaffolds
while read -r SCAFFOLD; do
    echo "Processing scaffold: $SCAFFOLD"
    
    VCF_IN="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.vcf.gz"
    VCF_OUT="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.phased.vcf.gz"

    if [[ ! -f "$VCF_IN" ]]; then
        echo "Warning: Input VCF file for scaffold $SCAFFOLD not found. Skipping."
        continue
    fi
    if [[ -f "$VCF_OUT" ]]; then
        echo "Phased VCF already exists for scaffold $SCAFFOLD; skipping."
        continue
    fi

    echo "Phasing VCF for scaffold $SCAFFOLD..."
    
    # Remove leading spaces from VCF
    #sed -i.bak 's/^ //g' "$VCF_IN"

    # Run shapeit phasing
    shapeit4 phase_common --input "$VCF_IN" --region 1 --output "$VCF_OUT" --thread $THREADS
    
    if [[ $? -ne 0 ]]; then
        echo "Error: shapeit phase failed for scaffold $SCAFFOLD."
        exit 1
    fi

    # Perform sapphire polishing

    if [ "$SAPPHIRE" = true ]; then

        if [ -d "${OUTDIR}/datafiles/sapphire_polishing" ]; then
        echo "Directory ${OUTDIR}/datafiles/sapphire_polishing exists."
        else
        echo "Directory ${OUTDIR}/datafiles/sapphire_polishing does not exist Creating it now."
        mkdir -p ${OUTDIR}/datafiles/sapphire_polishing
        fi

        echo "Polishing phased VCF for scaffold $SCAFFOLD with SAPPHIRE..."

        VCF_ANNOTATED="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.phased_PP_annotated.vcf"
        VCF_REPHASED="${OUTDIR}/datafiles/split_vcf/${POP}_${SCAFFOLD}.rephased.vcf"
        EXTRACTED_PP="${OUTDIR}/datafiles/sapphire_polishing/${POP}_${SCAFFOLD}.phased_PP_extract.bin"
        EXTRACTED_PP_ORIGINAL="${OUTDIR}/datafiles/sapphire_polishing/${POP}_${SCAFFOLD}.phased_PP_extract.bin.original"
        ANNOTATED_CSV="${OUTDIR}/datafiles/sapphire_polishing/${POP}_${SCAFFOLD}.phased_PP_annotated.samples.csv"
        PP_INFO="${OUTDIR}/datafiles/sapphire_polishing/pp_info.tsv"
        HEADER="${OUTDIR}/sapphire_polishing/header.txt"

        # Extract all variants with AF < 0.01 and replace heterozygous by "0.5" and homozygous by "."
        bcftools filter -i 'INFO/AF<0.01' $VCF_OUT -Ou | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' | sed -e 's/0|0/./g' -e 's/0|1/0.5/g' -e 's/1|0/0.5/g' -e 's/1|1/./g' > $PP_INFO
        # bgzip and tabix annotate the file
        bgzip $PP_INFO
        tabix -s1 -b2 -e2 $PP_INFO.gz

        if [ -f "${HEADER}" ];
            then
                echo "Header file already exists, moving on!"
            else
                echo '##FORMAT=<ID=PP,Number=1,Type=Float,Description="SHAPEIT Phasing Probability (PP) field">' > ${HEADER}
        fi


        bcftools annotate -a ${PP_INFO}.gz -h ${HEADER} -c CHROM,POS,ID,REF,ALT,FORMAT/PP \
        $VCF_OUT -Ob -o $VCF_ANNOTATED

        # Sanity check the annotated VCF
        bcftools view "$VCF_ANNOTATED" | less

        ${PROGDIR}/sapphire/pp_extractor/pp_extract -f $VCF_ANNOTATED --show-number -o $EXTRACTED_PP
        # Copy the file as SAPPHIRE will modify it, with this we will be able to compare
        cp "$EXTRACTED_PP" "$EXTRACTED_PP_ORIGINAL"

        bcftools query --list-samples $VCF_ANNOTATED | \
        awk '{print NR-1 "," $0 ",CRAM_PATH"}' > \
        $ANNOTATED_CSV

        ${PROGDIR}/sapphire/phase_caller/phase_caller \
        -f $VCF_ANNOTATED \
        -S $ANNOTATED_CSV \
        --cram-path-from-samples-file \
        -b "$EXTRACTED_PP" \
        -t $(nproc)

        #./sapphire/bin_tools/bin_diff -a "$EXTRACTED_PP_ORIGINAL" \
        #-b "$EXTRACTED_PP" \
        #-f $VCF_OUT \
        #-S $ANNOTATED_CSV \
        #--extra-info --more | less

        ${PROGDIR}/sapphire/pp_update/pp_update -f $VCF_ANNOTATED \
        -b $EXTRACTED_PP \
        -o $VCF_REPHASED


        
        #Remove the unpolished phased VCF and intermediate VCFs
        rm "$VCF_OUT"
        rm "$VCF_ANNOTATED"
        
    fi



done < "$SCAFFOLD_LIST"

echo "VCF phasing completed."
date
