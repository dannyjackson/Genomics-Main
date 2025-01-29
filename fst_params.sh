module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2
module load python/3.11/3.11.4

MSMCTOOLS=~/programs/msmc-tools # directory with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS

OUTDIR=/xdisk/mcnew/dannyjackson/cardinals_dfinch/    # main directory for output files
ANGSD=~/programs/angsd/     # path to directory with angsd executables

scriptdir=~/programs/Genomics-Main/
PATH=$PATH:$scriptdir # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path

# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/fst
mkdir -p ${OUTDIR}/analyses/genelist
mkdir -p ${OUTDIR}/datafiles/safs
mkdir -p ${OUTDIR}/datafiles/mls/
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
mkdir -p ${OUTDIR}/analyses/genelist/${WIN}


# for fst_1.sh

# define the names of the two populations that will be compared
POP1=nocaurban
POP2=nocarural

# define two colors to be used 
color1=#4EAFAF
color2=#082B64



# characters at the start of a chromosome number (excluding scaffolds)
# note that this script also assumes chromosomes will end in .1 i.e. NC_012345.1, and will remove the .1 from files where it will disrupt plotting etc

CHRLEAD=NC_0
SEXCHR=NC_044601
# path to reference genome
REF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
# path to gff file
GFF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCA_013397215.1/genomic.gff 



# make reference files
# first, make a file with chromosome name and length of chromosome
awk 'BEGIN {OFS = "\t"} {print $1,$2}' ${REF}.fai | grep ${CHRLEAD} ${OUTDIR}/referencelists/allscaffs_lengths.txt | grep -v ${SEXCHR} ${OUTDIR}/referencelists/allchroms_lengths.txt > ${OUTDIR}/referencelists/autosomes_lengths.txt

while IFS=',' read -r first second; do
    sed -i "s/$second/$first/g" ${OUTDIR}/referencelists/autosomes_lengths.txt 
done <<< "$CHROM"


# Make a comma separated chromosome conversion file without a header where the first column is the name of the chromosome and the second is the name of the associated scaffold in the reference genome:

if [ -f "${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt"]
        then
            echo "Chromosome conversion table already complete, moving on!"
        else
        echo '1,NC_044571.1' > ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '2,NC_044572.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '3,NC_044573.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '4,NC_044574.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '5,NC_044575.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '6,NC_044576.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '7,NC_044577.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '8,NC_044578.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '9,NC_044579.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '10,NC_044580.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '11,NC_044581.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '12,NC_044582.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '13,NC_044583.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '14,NC_044584.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '15,NC_044585.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '1A,NC_044586.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '17,NC_044587.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '18,NC_044588.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '19,NC_044589.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '20,NC_044590.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '21,NC_044591.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '22,NC_044592.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '23,NC_044593.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '24,NC_044594.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '25,NC_044595.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '26,NC_044596.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '27,NC_044597.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '28,NC_044598.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '29,NC_044599.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo '4A,NC_044600.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
        echo 'Z,NC_044601.1' >> ${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
fi