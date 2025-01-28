#!/bin/sh

# FST script (2/3)
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
    while getopts p:w:s: option
    do
    case "${option}"
    in
    p) PARAMS=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;

    esac
    done

source ${PARAMS}

WIN="${WIN:-10000}"
STEP="${STEP:-10000}"


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
grep ${CHRLEAD} ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2} >> ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
sed -i 's/${CHRLEAD}//g' ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
sed -i 's/\.1\t/\t/g' ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt

Rscript ${scriptdir}/fst_window.r ${OUTDIR} ${WIN} ${POP1} ${POP2}



if [ -f "${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}"* ]
        then
            echo "SNP analysis already complete, moving on!"
        else
            echo "computing fst on single snps"
            
            ${ANGSD}/misc/realSFS fst stats2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx -win 1 -step 1 >${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2}

            # single snps
            echo -e 'region\tchr\tmidPos\tNsites\tfst' > ${OUTDIR}/analyses/fst/singlesnps_fst_pyrr.txt
            #tail -n+2 slidingwindow >> slidingwindow_fst.txt 
            grep ${CHRLEAD} ${OUTDIR}/analyses/fst/singlesnps_pyrr >> ${OUTDIR}/analyses/fst/singlesnps_fst_pyrr.txt
            sed -i 's/${CHRLEAD}//g' ${OUTDIR}/analyses/fst/singlesnps_fst_pyrr.txt 
            sed -i 's/\.1\t/\t/g' ${OUTDIR}/analyses/fst/singlesnps_fst_pyrr.txt

            Rscript ${OUTDIR}/programs/Intro_Bioinformatics_Workshop/fst_snps.r ${OUTDIR} ${WIN}
fi


fi