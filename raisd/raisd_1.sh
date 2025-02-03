#!/bin/sh

# pausing on my development of these scripts, get msmc working first bc the input files are probably useful here

# raisd script (1/N)

if [ $# -lt 1 ]; then
    cat <<EOF
This script applies the software RAiSD to detect selective sweeps from whole-genome sequence data

It requires a bcf file as input, which can be generated using the <scriptname> scripts in:
    github.com/dannyjackson/Genomics-Main

Read and understand the entire script before running it!

REQUIRED ARGUMENTS:
    -p  Path to parameter file (example: params_raisd.sh in GitHub repo)
EOF
    exit 1
fi

# Parse arguments
while getopts "p:w:s:c:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        s) STEP=${OPTARG} ;;
        c) CHR_FILE=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

source "${PARAMS}"

if [ -z "${OUTDIR}" ] || [ -z "${POP1}" ] || [ -z "${POP2}" ]; then
    echo "Error: OUTDIR, POP1, or POP2 is not set in the parameter file." >&2
    exit 1
fi

if [ ! -f "${CHR_FILE}" ]; then
    echo "Error: Chromosome mapping file ${CHR_FILE} not found." >&2
    exit 1
fi

CHROM=`cat $CHR_FILE`

printf "\n \n \n \n"
date
echo "Current script: dxy.sh"

cd "${OUTDIR}/analyses/dxy/${POP1}_${POP2}" || { echo "Failed to enter directory"; exit 1; }

# B5_raisd
# for setup file
mkdir -p ${OUTDIR}/analyses/raisd
mkdir -p ${OUTDIR}/analyses/raisd/${POP1}
mkdir -p ${OUTDIR}/analyses/raisd/${POP2}

# subset vcf by autosomes
# subset vcf by population
# convert to bcf

bcftools convert /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/noca_urban.bcf \
    -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/noca_urban.vcf -O v 


# install raisd
mkdir RAiSD
cd RAiSD
wget https://github.com/alachins/raisd/archive/master.zip
unzip master.zip
cd raisd-master
./install-RAiSD.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/gsl/lib


/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n /xdisk/mcnew/dannyjackson/cardinals/raisd/output/noca_urban -I /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/noca_urban.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna


