#!/bin/sh 

# raisd script

# RAiSD script for detecting selective sweeps from whole-genome sequence data.

# Function to display usage information
usage() {
    cat <<EOF
This script applies the software RAiSD to detect selective sweeps from whole-genome sequence data.

It requires a gzipped vcf file as input, which can be generated using the <scriptname> scripts in:
    github.com/dannyjackson/Genomics-Main

REQUIRED ARGUMENTS:
    -p  path to params file
    -n  population name
    -w  Window size
EOF
    exit 1
}


# Parse command-line arguments
while getopts "p:n:w:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        n) POP=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

source "${PARAMS}"

cd ${OUTDIR}/analyses/raisd/${POP}/${WIN}

while read -r SCAFFOLD; do
    VCF_IN=/xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/${POP}.${SCAFFOLD}.phased.vcf

    ~/programs/RAiSD/raisd-master/RAiSD \
        -n "${POP}.${SCAFFOLD}" \
        -I "${VCF_IN}" \
        -f -O -R -P -C "${REF}" -w ${WIN}
        
done < "$SCAFFOLD_LIST"

mv /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/RAiSD_Info* /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/infofiles

mv /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/RAiSD_Plot* /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/plots

mv /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/RAiSD_Report* /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${POP}/${WIN}/reportfiles

echo -e "chrom\tposition\tstart\tend\tVAR\tSFS\tLD\tU" > ${OUTDIR}/analyses/raisd/${POP}/${WIN}/RAiSD_Report.${POP}.chromosomes

while read -r SCAFFOLD;
do 
  awk -v CHROM="$SCAFFOLD" '{print CHROM, $0}' ${OUTDIR}/analyses/raisd/${POP}/${WIN}/reportfiles/RAiSD_Report.${POP}.${SCAFFOLD}.${SCAFFOLD} >> ${OUTDIR}/analyses/raisd/${POP}/${WIN}/RAiSD_Report.${POP}.chromosomes

done < "${OUTDIR}/referencelists/SCAFFOLDS.txt"


# rename scaffolds for plotting

echo "Processing CHROM file: $CHR_FILE..."

FILE="${OUTDIR}/analyses/raisd/${POP}/${WIN}/RAiSD_Report.${POP}.chromosomes"
# Read CHROM line by line
while IFS=',' read -r first second; do
    echo "Replacing occurrences of '$second' with '$first' in $FILE"
    sed -i.bak "s/$second/$first/g" "$FILE"
done < "$CHR_FILE"

rm -f "${FILE}.bak"

# z transform U metric
WIN_OUT="${OUTDIR}/analyses/${METRIC}/${POP}/${POP}.${METRIC}_${WIN}.Ztransformed.csv"

echo 'visualizing windows'
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}"
echo 'finished windowed plot' 

# plot scaffolds
echo 'visualizing windows'
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "${POP}"
echo 'finished windowed plot' 