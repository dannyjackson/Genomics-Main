module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2
module load python/3.11/3.11.4
module load bwa/0.7.17
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module load samtools/1.19.2

# Define variables
# all
OUTDIR=/xdisk/mcnew/dannyjackson/cardinals_dfinch/    # main directory for output files
PROGDIR=~/programs
BAMDIR=/path/to/bam/files
PROJHUB=github_project_name
scriptdir=${PROGDIR}/${PROJHUB}
PATH=$PATH:$scriptdir # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path
THREADS=4
ID=name_of_project
FILENAME_LIST="/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt" # list with sample codes associated with each file in dataset, one per line

# define the names of the two populations that will be compared
POP1=nocaurban
POP2=nocarural
# define two colors to be used 
color1=#4EAFAF
color2=#082B64
# note that this script also assumes chromosomes will end in .1 i.e. NC_012345.1, and may also remove the .1 from files where it will disrupt plotting etc
CHRLEAD=NC_0 # characters at the start of a chromosome number (excluding scaffolds)
SEXCHR=NC_044601
REF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna # path to reference genome
GFF=/xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCA_013397215.1/genomic.gff # path to gff file

# Generate scaffold list
if [ -f "${scriptdir}/SCAFFOLDS.txt"]
        then
            echo "SCAFFOLDS.txt already exists, moving on!"
        else
        awk '{print $1}' "${GENOME}.fai" > "${scriptdir}/SCAFFOLDS.all.txt"
        grep "$CHRLEAD" "${scriptdir}/SCAFFOLDS.all.txt" > "${scriptdir}/SCAFFOLDS.chroms.txt"
        grep -v "$SEXCHR" "${scriptdir}/SCAFFOLDS.chroms.txt" > "${scriptdir}/SCAFFOLDS.txt"
fi

# specific to selection analyses (fst, dxy, Tajima's D, RAiSD)
ANGSD=~/programs/angsd/ # path to directory with angsd executables

# msmc specific
k=150 # For snpability.sh
prefix=GCF_901933205
MSMCTOOLS=${PROGDIR}/msmc-tools # directory with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS # add directory with msmc-tools binaries to path
METHOD=samtools  # or another variant calling method


# other code
# specific to selection analyses (fst, dxy, Tajima's D, RAiSD)
# make directories for intermediate files-- will fail if these don't exist

mkdir -p ${OUTDIR}/analyses/fst
mkdir -p ${OUTDIR}/analyses/genelist
mkdir -p ${OUTDIR}/datafiles/safs
mkdir -p ${OUTDIR}/datafiles/mls/
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
mkdir -p ${OUTDIR}/analyses/genelist/${WIN}












# Miscellaneous
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


