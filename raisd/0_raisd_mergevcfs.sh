#!/bin/bash

module load bcftools
module load samtools

for GROUP in nocaurban nocarural pyrrurban pyrrrural; do
    while read -r IND; do
        while read -r SCAFFOLD; do
            echo "Indexing VCF for: $IND $SCAFFOLD"
            bcftools index /xdisk/mcnew/dannyjackson/cardinals/datafiles/vcf2/${IND}.${SCAFFOLD}.phased.vcf.gz
        done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt
    done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_${GROUP}.txt
done

for GROUP in nocaurban nocarural pyrrurban pyrrrural; do

    while read -r SCAFFOLD; do
        echo "Merging VCFs for scaffold: $SCAFFOLD"
        
        # Create a list of input VCF files
        INPUT_VCFS=()
        while read -r IND; do
            FILE="/xdisk/mcnew/dannyjackson/cardinals/datafiles/vcf2/${IND}.${SCAFFOLD}.phased.vcf.gz"
            if [[ -f "$FILE" ]]; then
                INPUT_VCFS+=("$FILE")
            else
                echo "Warning: File $FILE not found."
            fi
        done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_${GROUP}.txt

        # Merge if there are files to merge
        if [[ ${#INPUT_VCFS[@]} -gt 0 ]]; then
            bcftools merge "${INPUT_VCFS[@]}" -O v -o /xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/${GROUP}.${SCAFFOLD}.phased.vcf
        else
            echo "No valid VCF files found for scaffold $SCAFFOLD."
        fi
    done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt
done 
