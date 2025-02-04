#!/bin/sh

# tajima script 

if [ $# -lt 1 ]; then
    cat <<EOF
This script computes Tajima's D within a population of genomes using a genotype likelihood framework implemented in ANGSD.

It requires SAF files as input, which can be generated using the <scriptname> scripts in:
    github.com/dannyjackson/Genomics-Main

It computes average genome-wide tajima's d and produces output files for:
    - Sliding window dxy
    - Per-SNP dxy

Read and understand the entire script before running it!

REQUIRED ARGUMENTS:
    -p  Path to parameter file (example: params_dxy.sh in GitHub repo)
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

if [ -z "${OUTDIR}" ] || [ -z "${POP}" ] || [ -z "${WIN}"] || [ -z "${STEP}" ]; then
    echo "Error: OUTDIR, POP, WIN, or STEP is not set in the parameter file." >&2
    exit 1
fi


printf "\n \n \n \n"
date
echo "Current script: tajima.sh"



echo "Generating sfs file"
${ANGSD}/misc/realSFS -P 24 ${OUTDIR}/datafiles/safs/${POP}.saf.idx > ${OUTDIR}/datafiles/safs/${POP}.sfs 

echo "Calculating thetas for each site"
${ANGSD}/misc/realSFS saf2theta ${OUTDIR}/datafiles/safs/${POP}.saf.idx -outname ${OUTDIR}/analyses/thetas/${POP} -sfs ${OUTDIR}/datafiles/safs/${POP}.sfs

echo "Estimating Tajima's D genome wide"
${ANGSD}/misc/thetaStat do_stat ${OUTDIR}/analyses/thetas/${POP}.thetas.idx

echo "Estimating sliding window Tajima's D"
${ANGSD}/misc/thetaStat do_stat ${OUTDIR}/analyses/thetas/${POP}.thetas.idx -win ${WIN} -step ${STEP}  -outnames ${OUTDIR}/analyses/thetas/${POP}/${WIN}/${POP}.theta.thetasWindow

WIN_OUT="${OUTDIR}/analyses/thetas/${POP}/${WIN}/${POP}.theta.thetasWindow.pestPG"

# Check if CHROM has anything assigned
if [[ -n "$CHROM" ]]; then
    echo "Processing CHROM variable..."


    # Read CHROM line by line
    while IFS=',' read -r first second; do
        echo "Replacing occurrences of '$second' with '$first' in $WIN_OUT"
        sed -i.bak "s/$second/$first/g" "$WIN_OUT"
    done <<< "$CHROM"

    rm -f "${WIN_OUT}.bak"
else
    echo "CHROM variable is empty or not set."
fi

Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}" 

