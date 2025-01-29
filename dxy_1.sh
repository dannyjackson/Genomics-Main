#!/bin/sh

# dxy script (1/N)

if [ $# -lt 1 ]
  then
    echo "This script computes dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide dxy and produce the output files necessary for sliding window dxy and dxy for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as fst_params.sh"

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

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"

CHROM=`cat $CHR_FILE`

printf "\n \n \n \n"
date
echo "Current script: dxy_1.sh"

cd ${OUTDIR}/analyses/dxy/${POP1}_${POP2}

# unzip maps files if necessary
# Define file paths
POP1_MAFS_GZ="${OUTDIR}/datafiles/safs/${POP1}.mafs.gz"
POP1_MAFS="${OUTDIR}/datafiles/safs/${POP1}.mafs"

# Check if the .gz file exists
if [ -f "$POP1_MAFS_GZ" ]; then
    # If .gz file exists, unzip it
    echo "Found ${POP1}.mafs.gz, unzipping..."
    gunzip "$POP1_MAFS_GZ"
elif [ -f "$POP1_MAFS" ]; then
    # If .gz file doesn't exist, but the .mafs file exists, continue
    echo "Found ${POP1}.mafs, continuing..."
else
    # If neither file exists, break
    echo "Neither ${POP1}.mafs.gz nor ${POP1}.mafs found. Exiting..."
    exit 1
fi

# Define file paths
POP2_MAFS_GZ="${OUTDIR}/datafiles/safs/${POP2}.mafs.gz"
POP2_MAFS="${OUTDIR}/datafiles/safs/${POP2}.mafs"

# Check if the .gz file exists
if [ -f "$POP2_MAFS_GZ" ]; then
    # If .gz file exists, unzip it
    echo "Found ${POP2}.mafs.gz, unzipping..."
    gunzip "$POP2_MAFS_GZ"
elif [ -f "$POP2_MAFS" ]; then
    # If .gz file doesn't exist, but the .mafs file exists, continue
    echo "Found ${POP2}.mafs, continuing..."
else
    # If neither file exists, break
    echo "Neither ${POP2}.mafs.gz nor ${POP2}.mafs found. Exiting..."
    exit 1
fi


total_lines=$(cat ${OUTDIR}/datafiles/safs/${POP1}.mafs | wc -l)
num_sites=$((total_lines - 1))

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p ${OUTDIR}/datafiles/safs/${POP1}.mafs -q ${OUTDIR}/datafiles/safs/${POP2}.mafs -t ${num_sites} > ${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_globalestimate_${POP1}_${POP2}.txt

mv ${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite.txt ${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.txt

# Check if CHROM has anything assigned
if [[ -n "$CHROM" ]]; then
    echo "Processing CHROM variable..."
    
    # Define the files to process
    files=(
        "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.txt"
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


awk '{print $2}' ${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.txt  | tail -n +2 > ${OUTDIR}/analyses/dxy/${POP1}_${POP2}/${POP1}_${POP2}_sites.txt


Rscript ~/programs/Genomics-Main/dxy_snps.r ${OUTDIR} ${POP1} ${POP2} ${COLOR1} ${COLOR2}

python ~/programs/Genomics-Main/dxy_windows.py --outdir ${OUTDIR} --pop1 ${POP1} --pop2 ${POP2} --win ${WIN}

fi