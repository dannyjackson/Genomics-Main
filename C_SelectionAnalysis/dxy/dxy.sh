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


printf "\n \n \n \n"
date

echo "Current script: dxy.sh"

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
# Define output file paths
DXY_OUTPUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_globalestimate_${POP1}_${POP2}.txt"
PERSITE_OUTPUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/Dxy_persite_${POP1}_${POP2}.txt"
AUTOSOMES_OUTPUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/Dxy_persite_${POP1}_${POP2}.autosomes.txt"
SITES_OUTPUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/${POP1}_${POP2}_sites.txt"

# Run Dxy calculation only if the output file does not already exist
if [[ ! -f "$DXY_OUTPUT" ]]; then
    Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R \
        -p "${OUTDIR}/datafiles/safs/${POP1}.mafs" \
        -q "${OUTDIR}/datafiles/safs/${POP2}.mafs" \
        -t "${num_sites}" > "$DXY_OUTPUT"
else
    echo "Dxy calculation already completed. Skipping."
fi

# Move the per-site Dxy file only if it doesn't exist
if [[ -f "$PERSITE_OUTPUT" ]]; then
    mv "$PERSITE_OUTPUT" \
       "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/Dxy_persite_${POP1}_${POP2}.txt"
else
    echo "Per-site Dxy file not found or already moved."
fi

# Write header and filter autosomes only if the output file doesn't already exist
if [[ ! -f "$AUTOSOMES_OUTPUT" ]]; then
    echo -e "chromo\tposition\tdxy" > "$AUTOSOMES_OUTPUT"
    grep ${CHRLEAD} "$PERSITE_OUTPUT" | grep -v ${SEXCHR} >> "$AUTOSOMES_OUTPUT"
else
    echo "Autosomes filtered file already exists. Skipping."
fi

# Process CHROM file only if it exists
if [[ -f "$CHR_FILE" ]]; then
    echo "Processing CHROM file: $CHR_FILE..."

    FILE="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}_average_dxy_${WIN}bp_windows.txt"

    if [[ -f "$FILE" ]]; then
        while IFS=',' read -r first second; do
            echo "Replacing occurrences of '$second' with '$first' in $FILE"
            sed -i.bak "s/$second/$first/g" "$FILE"
        done < "$CHR_FILE"

        rm -f "${FILE}.bak"
    else
        echo "Windowed Dxy file not found. Skipping CHROM processing."
    fi
else
    echo "CHROM file not set or missing. Skipping CHROM processing."
fi

# Extract site positions only if the output file doesn't already exist
if [[ ! -f "$SITES_OUTPUT" ]]; then
    awk 'NR>1 {print $2}' "$AUTOSOMES_OUTPUT" > "$SITES_OUTPUT"
else
    echo "Site positions file already exists. Skipping."
fi


# Compute windows and produce manhattan plots for windows and snp data
# Define output files
WIN_OUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}_average_dxy_${WIN}bp_windows.txt"
SNP_IN="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/Dxy_persite_${POP1}_${POP2}.autosomes.txt"
SNP_OUT="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/snps/${POP1}_${POP2}.dxy.snps.sigline.png"
       
# Run first two scripts in sequence if output file doesn't exist
if [ ! -f "$WIN_OUT" ]; then
    echo 'computing windows'
    python "${SCRIPTDIR}/Genomics-Main/dxy/dxy_windows.py" --outdir "${OUTDIR}" --pop1 "${POP1}" --pop2 "${POP2}" --win "${WIN}" 

    # Check if CHROM has anything assigned
if [[ -f "$CHR_FILE" ]]; then
    echo "Processing CHROM file: $CHR_FILE..."
    
    # Define the files to process

    FILE="${OUTDIR}/analyses/dxy/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}_average_dxy_${WIN}bp_windows.txt"


    # Read CHROM line by line
    while IFS=',' read -r first second; do
        echo "Replacing occurrences of '$second' with '$first' in $FILE"
        sed -i.bak "s/$second/$first/g" "$FILE"
    done < "$CHR_FILE"

    rm -f "${FILE}.bak"
else
    echo "CHROM variable is empty or not set."
fi
    echo 'visualizing windows'
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP1}" "${POP2}"
    echo 'finished windowed plot' &
fi

# Run the SNP visualization separately if output file doesn't exist
if [ ! -f "$SNP_OUT" ]; then
    echo 'visualizing snps'
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${SNP_IN}" "snps" "${POP1}" "${POP2}"
    echo 'finished snp plot' &
fi

# Wait for background jobs to finish
wait


#Rscript ${SCRIPTDIR}/Genomics-Main/dxy/dxy_snps.r "${OUTDIR}" "${POP1}" "${POP2}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${OUTDIR}/analyses/dxy/${POP1}_${POP2}/Dxy_persite_nocaurban_nocarural.txt"




# Rscript ${SCRIPTDIR}/Genomics-Main/dxy/dxy_windows.r "${OUTDIR}" "${POP1}" "${POP2}" "${COLOR1}" "${COLOR2}" "${CUTOFF}"


