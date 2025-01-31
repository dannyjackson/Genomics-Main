#!/bin/sh

# FST script (3/4)
# This file accesses fst_window.r and fst_snps.r

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst and produce the output files necessary for sliding window Fst and fst for each SNP. 
    
    Read and understand the entire script before running it! 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as params.sh)

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)"

  else
    while getopts p:w:s:c: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;
    c) CHR_FILE=${OPTARG};;

    esac
    done

source ${PARAMS}

WIN="${WIN:-10000}"
STEP="${STEP:-10000}"
CHROM=`cat $CHR_FILE`

printf "\n \n \n \n"
date
echo "Current script: fst_2.sh"


# sliding window analyses
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
# windowed
if [ -f "${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}"* ]
        then
            echo "windowed fst output file already present, assuming it is already generated and moving on!"
        else
            # sliding window analysis
            echo "computing sliding window statistics"
            ${ANGSD}/misc/realSFS fst stats2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx -win ${WIN} -step ${STEP} > ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}

fi
 


# windowed
echo -e 'region\tchr\tmidPos\tNsites\tfst' > ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
grep ${CHRLEAD} ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2} | grep -v ${SEXCHR} >> ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt


# replace chromosome names if necessary

# CHROM="/xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt"

# Check if CHROM has anything assigned
if [[ -n "$CHROM" ]]; then
    echo "Processing CHROM variable..."
    
    # Define the files to process
    files=(
        "$OUTDIR/analyses/fst/$WIN/slidingwindow.${POP1}_${POP2}"
        "$OUTDIR/analyses/fst/$WIN/slidingwindow.${POP1}_${POP2}.chroms.txt"
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

fi
Rscript ${scriptdir}/Genomics-Main/fst_window.r ${OUTDIR} ${WIN} ${POP1} ${POP2}
