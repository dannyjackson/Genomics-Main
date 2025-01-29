#!/bin/sh

# FST script (2/4)
# This file accesses fst_window.r and fst_snps.r

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst and produce the output files necessary for sliding window Fst and fst for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as fst_params.sh)

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)"

  else
    while getopts p:c: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;
    c) CHR_FILE=${OPTARG};;

    esac
    done

source ${PARAMS}

CHROM=`cat $CHR_FILE`

printf "\n \n \n \n"
date
echo "Current script: fst_2.sh"

if [ -f "${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}"* ]
        then
            echo "SNP analysis already complete, moving on!"
        else
            echo "computing fst on single snps"
            
            ${ANGSD}/misc/realSFS fst stats2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx -win 1 -step 1 >${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}

            # single snps
            echo -e 'region\tchr\tmidPos\tNsites\tfst' > ${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}.chroms.txt 
            #tail -n+2 slidingwindow >> slidingwindow_fst.txt 
            grep ${CHRLEAD} ${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2} | grep -v ${SEXCHR} >> ${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}.chroms.txt

            sed -i 's/${CHRLEAD}//g' ${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}.chroms.txt 

            # Check if CHROM has anything assigned
            if [[ -n "$CHROM" ]]; then
                echo "Processing CHROM variable..."
                
                # Define the files to process
                files=(
                    "$OUTDIR/analyses/fst/singlesnps.${POP1}_${POP2}"
                    "$OUTDIR/analyses/fst/singlesnps.${POP1}_${POP2}.chroms.txt"
                )

                # Read CHROM line by line
                while IFS=',' read -r first second; do
                    echo "Replacing occurrences of '$second' with '$first'..."
                    
                    # Process each file
                    for file in "${files[@]}"; do
                        if [[ -f "$file" ]]; then
                            echo "Processing file: $file"
                            sed -i "s/$second/$first/g" "$file"
                        else
                            echo "Warning: File $file not found."
                        fi
                    done
                done <<< "$CHROM"

            else
                echo "CHROM variable is empty or not set."
            fi
            Rscript ${scriptdir}/fst_snps.r ${OUTDIR} ${POP1} ${POP2} ${COLOR1} ${COLOR2}
fi
fi