#!/bin/sh

# dxy script (1/N)

if [ $# -lt 1 ]; then
    cat <<EOF
This script computes dxy between two groups of genomes using a genotype likelihood framework implemented in ANGSD.
It requires SAF files as input, which can be generated using the <scriptname> scripts in:
    github.com/dannyjackson/Genomics-Main

It computes average genome-wide dxy and produces output files for:
    - Sliding window dxy
    - Per-SNP dxy

**NOTE:** Read the entire script and adjust for each project! Many parameters (e.g., SNP filtering in ANGSD) are hardcoded.

REQUIRED ARGUMENTS:
    -p  Path to parameter file (example: params.sh in GitHub repo)
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
echo "Current script: dxy_1.sh"

cd "${OUTDIR}/analyses/dxy/${POP1}_${POP2}" || { echo "Failed to enter directory"; exit 1; }

# unzip maps files if necessary
for POP in "${POP1}" "${POP2}"; do
    MAF_GZ="${OUTDIR}/datafiles/safs/${POP}.mafs.gz"
    MAF="${OUTDIR}/datafiles/safs/${POP}.mafs"

    if [ -f "${MAF_GZ}" ]; then
        echo "Found ${POP}.mafs.gz, unzipping..."
        gzip -d -f "${MAF_GZ}"
    elif [ -f "${MAF}" ]; then
        echo "Found ${POP}.mafs, continuing..."
    else
        echo "Error: Neither ${POP}.mafs.gz nor ${POP}.mafs found. Exiting..."
        exit 1
    fi
done

# Calculate number of sites
total_lines=$(wc -l < "${OUTDIR}/datafiles/safs/${POP1}.mafs")
num_sites=$((total_lines - 1))

# Run Dxy calculation
Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R \
    -p "${OUTDIR}/datafiles/safs/${POP1}.mafs" \
    -q "${OUTDIR}/datafiles/safs/${POP2}.mafs" \
    -t "${num_sites}" > "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_globalestimate_${POP1}_${POP2}.txt"

mv "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite.txt" \
   "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.txt"

# Write header to the output file
echo -e "chromo\tposition\tdxy" > "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.autosomes.txt"

# Filter the input file, excluding sex chromosomes, and append results
grep ${CHRLEAD} "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.txt" | grep -v ${SEXCHR} >> "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.autosomes.txt"

# Check if CHROM has anything assigned
if [[ -n "$CHROM" ]]; then
    echo "Processing CHROM variable..."
    
    # Define the files to process
    FILE="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.autosomes.txt"

    # Read CHROM line by line
    while IFS=',' read -r first second; do
        echo "Replacing occurrences of '$second' with '$first' in $FILE"
        sed -i.bak "s/$second/$first/g" "$FILE"
    done <<< "$CHROM"

    rm -f "${FILE}.bak"
else
    echo "CHROM variable is empty or not set."
fi

# Extract site positions
awk 'NR>1 {print $2}' "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_${POP1}_${POP2}.autosomes.txt" \
    > "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/${POP1}_${POP2}_sites.txt"


# Run visualization scripts
Rscript ${PROGDIR}/Genomics-Main/dxy_snps.r "${OUTDIR}" "${POP1}" "${POP2}" "${COLOR1}" "${COLOR2}"

python ${PROGDIR}//Genomics-Main/dxy_windows.py --outdir "${OUTDIR}" --pop1 "${POP1}" --pop2 "${POP2}" --win "${WIN}"