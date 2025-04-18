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
    -w  window size
    -r reference genome
EOF
    exit 1
}


# Parse command-line arguments
while getopts "p:n:w:" option; do
    case "${option}" in
        p) PARAMS=${OPTARG} ;;
        n) POP=${OPTARG} ;;
        w) WIN=${OPTARG} ;;
        r) REF=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

if [ -z "${PARAMS}" ]; then
    echo "Error: No parameter file provided." >&2
    exit 1
fi

source "${PARAMS}"

cd ${OUTDIR}/analyses/raisd/${POP}/${WIN}


VCF_IN=${OUTDIR}/datafiles/vcf3/${POP}.phased.vcf

~/programs/RAiSD/raisd-master/RAiSD \
    -n "${POP}" \
    -I "${VCF_IN}" \
    -f -O -R -P -C "${REF}" -w ${WIN}
    

