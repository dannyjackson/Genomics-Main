module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2
module load python/3.11/3.11.4
module load bwa/0.7.18
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module load samtools/1.19.2



# Define variables
# all
OUTDIR=/path/to/project/directory/    # main directory for output files
PROGDIR=~/programs  # path to directory for all installed programs
BAMDIR=/path/to/bam/files
PROJHUB=github_project_name
SCRIPTDIR=${PROGDIR}/${PROJHUB}
PATH=$PATH:$SCRIPTDIR # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path
ID=name_of_project
FILENAME_LIST="/path/to/list.txt" # list with sample codes associated with each file in dataset, one per line

# define aspects of the reference genome
CHRLEAD=NC_0 # characters at the start of a chromosome number (excluding scaffolds)
SEXCHR=NC_044601
REF=/path/to/reference/genome/file.fna # path to reference genome
GFF=/path/to/reference/genome/gff/genomic.gff # path to gff file
