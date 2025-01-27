#!/bin/sh

# FST script
# This file accesses fst_window.r and fst_snps.r

if [ $# -lt 1 ]
  then
    echo "This script computes Fst between two groups of genomes using a genotype likelihood framework implemented in ANGSD. It requires SAF files as input, which can be generated using the <scriptname> scripts in github.com/dannyjackson/Genomics-Main. It will compute average genome-wide Fst, sliding window Fst, and fst for each SNP. Read the entire script and revise according to each project! Many parameters are not modifiable with options, including snp filtering settings in angsd. 

    If you have run it once, feel free to change the window size and rerun. It will skip redundant analyses if the output files already exist (including Fst at individual SNPs) and only compute windowed analyses.

    [-o] Output directory
    [-p] Path to scripts
    [-x] Name of population 1
    [-y] Name of population 2
    [-r] Path to reference genome fasta
    [-g] Path to reference genome gff file

    OPTIONAL ARGUMENTS

    [-w] Window size for Fst scans (defaults to 10,000)
    [-s] Step size for Fst scans (defaults to 10,000)
    [-a] Path to angsd directory, will install angsd here if it is not already installed! (defaults to ~/program/angsd/)
    [-c] Leading text before numbers of chromosomes (include only chromosomes not scaffolds; defaults to "NC_0")

  else
    while getopts i:r:t:p:b:s: option
    do
    case "${option}"
    in
    o) OUTDIR=${OPTARG};;
    p) scriptdir=${OPTARG};;
    x) POP1=${OPTARG};;
    y) POP2=${OPTARG};;
    w) WIN=${OPTARG};;
    s) STEP=${OPTARG};;
    a) ANGSD=${OPTARG};;
    c) CHRLEAD=${OPTARG};;
    esac
    done

WIN="${WIN:-10000}"
STEP="${STEP:-10000}"
ANGSD="${ANGSD:-~/programs/angsd/}"
CHRLEAD="${CHRLEAD:-NC_0}"

module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2

PATH=$PATH:$scriptdir # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path

# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/fst
mkdir -p ${OUTDIR}/analyses/genelist
mkdir -p ${OUTDIR}/datafiles/safs
mkdir -p ${OUTDIR}/datafiles/mls/
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
mkdir -p ${OUTDIR}/analyses/genelist/${WIN}




printf "\n \n \n \n"
date
echo "Current script: fst.sh"

if [ -f "${ANGSD}" ]
        then
                echo "angsd already installed, moving on!"
        else
                echo "Install angsd"
                cd ~/IntroBioinformaticsWorkshop/programs/
                git clone https://github.com/ANGSD/angsd.git 
                cd angsd 
                make 
fi

# Generate saf files for each population in angsd. Skip if there is any output already in ${OUTDIR}/datafiles/safs/.

if [ -f "${OUTDIR}/datafiles/safs/" ]
        then
            echo "Files present in safs directory, assuming they are already generated and moving on!"
        else
            echo "computing safs for pop 1"
            ~/programs/angsd/angsd -bam ${OUTDIR}/datafiles/${POP1}.bams.txt -out ${OUTDIR}/datafiles/safs/${POP1} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites ${OUTDIR}/referencelists/sites_headless.mafs -anc ${REF} -nThreads 10
            
            echo "computing safs for pop 2"
            ~/programs/angsd/angsd -bam ${OUTDIR}/datafiles/${POP2}.bams.txt -out ${OUTDIR}/datafiles/safs/${POP1} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites ${OUTDIR}/referencelists/sites_headless.mafs -anc ${REF} -nThreads 10
fi


if [ -f "${OUTDIR}/datafiles/mls/" ]
        then
            echo "Files present in mls directory, assuming they are already generated and moving on!"
        else
            echo "computing the 2dsfs prior"
            ${ANGSD}/misc/realSFS ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx > ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml

fi


if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}*" ]
        then
            echo "fst.idx file already present, assuming it is already generated and moving on!"
        else
            echo "prepare the fst for easy window analysis, compute fst.idx file"
            ${ANGSD}/misc/realSFS fst index ${OUTDIR}/datafiles/safs/${POP1}.saf.idx ${OUTDIR}/datafiles/safs/${POP2}.saf.idx -sfs ${OUTDIR}/datafiles/mls/${POP1}_${POP2}.ml -fstout ${OUTDIR}/analyses/fst/${POP1}_${POP2}
fi

if [ -f "${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt" ]
        then
            echo "global fst estimate file already present, assuming it is already generated and moving on!"
        else
            echo "computing global fst estimate"
            ${ANGSD}/misc/realSFS fst stats ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt

fi

# sliding window analyses
# windowed
if [ -f "${OUTDIR}/analyses/fst/${WIN}/slidingwindow_${POP1}_${POP2}*" ]
        then
            echo "windowed fst output file already present, assuming it is already generated and moving on!"
        else
            # sliding window analysis
            echo "computing sliding window statistics"
            ${ANGSD}/misc/realSFS fst stats2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx -win ${WIN} -step ${STEP} > ${OUTDIR}/analyses/fst/${WIN}/slidingwindow_${POP1}_${POP2}

fi



# windowed
echo -e 'region\tchr\tmidPos\tNsites\tfst' > ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
grep ${CHRLEAD} ${OUTDIR}/analyses/fst/${WIN}/${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2} >> ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
sed -i 's/${CHRLEAD}//g' ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt
sed -i 's/\.1\t/\t/g' ${OUTDIR}/analyses/fst/${WIN}/slidingwindow.${POP1}_${POP2}.chroms.txt

Rscript ${scriptdir}/fst_window.r ${OUTDIR} ${WIN} ${POP1} ${POP2}



if [ -f "${OUTDIR}/analyses/fst/singlesnps.${POP1}_${POP2)*" ]
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



# make relevant gene lists


echo "Making relevant gene lists"

# snps 

awk 'BEGIN {FS = ","} {$1=""}1' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv.tmp

mv ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv.tmp ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv

sed -i 's/\"//g' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv

tail -n +2 ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.csv  > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.headless.csv


awk -F',' 'NR>1 {print $CHRLEAD, $1, ".1" "\t" $2-1 "\t" $2}' ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.headless.csv > ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.bed

bedtools intersect -a ${GFF} -b ${OUTDIR}/analyses/fst/${POP1}_${POP2}.outlierfst.bed -wa > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenes_snps_top.95.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenes_snps_top.95.txt | sed 's/ID\=gene\-//g' | sort -u > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.fst.relevantgenenames_snps_top.95.txt

# windowed

awk 'BEGIN {FS = ","} {$1=""}1' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv.tmp

mv ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv.tmp ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv

sed -i 's/\"//g' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv

tail -n +2 ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.csv > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.headless.csv


awk -F',' -v win="$WIN" 'NR>1 {print $CHRLEAD, $1, ".1" "\t" $2-(win/2) "\t" $2+(win/2)}' ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.windowed.outlierfst.headless.csv > ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.outlierfst.bed

bedtools intersect -a ${GFF} -b ${OUTDIR}/analyses/fst/${WIN}/${POP1}_${POP2}.outlierfst.bed -wa > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenes_windowed_top.95.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenes_windowed_top.95.txt | sed 's/ID\=gene\-//g' | sort -u > ${OUTDIR}/analyses/genelist/${POP1}_${POP2}.relevantgenenames_windowed_top.95.txt
