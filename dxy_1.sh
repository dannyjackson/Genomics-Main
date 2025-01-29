#!/bin/sh

# dxy script (1/N)

if [ $# -lt 1 ]
  then
    echo "This script computes dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide dxy and produce the output files necessary for sliding window dxy and dxy for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    REQUIRED ARGUMENTS
    [-p] Path to parameter file (example is saved in github repository as fst_params.sh"

  else
    while getopts p:w:s: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;

    esac
    done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi
source "${PARAMS}"


printf "\n \n \n \n"
date
echo "Current script: dxy_1.sh"


mkdir -p ${OUTDIR}/analyses/dxy/
mkdir -p ${OUTDIR}/analyses/dxy/${POP1}_${POP2}

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

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p ${OUTDIR}/datafiles/safs/${POP1}.mafs -q ${OUTDIR}/datafiles/safs/${POP2}.mafs -t ${num_sites}

mv Dxy_persite.txt Dxy_persite_${POP1}_${POP2}.txt

awk '{print $2}' Dxy_persite_${POP1}_${POP2}.txt  | tail -n +2 > ${POP1}_${POP2}_sites.txt

tail -n +2 Dxy_persite_${POP1}_${POP2}.txt > Dxy_persite_${POP1}_${POP2}.forplot.txt

Rscript ${scriptdir}/plotDXY.R -i Dxy_persite_${POP1}_${POP2}.forplot.txt -o ${POP1}_${POP2}_windowed -p ${POP1}_${POP2}_sites.txt -w ${WIN} -s ${STEP}

fi