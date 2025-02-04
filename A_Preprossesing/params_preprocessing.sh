module load bwa/0.7.17
module load samtools/1.19.2
module load bowtie2
module load picard
module load samtools
module load parallel
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9

source /path/to/params_base.sh
source /path/to/preprocessing_setup.sh

# across all preprocessing
THREADS=<SET_VALUE>

# trimming
FASTAS=/path/to/fasta/files # format must be samplename_
TRIMJAR=/path/to/trimmomatic/jarfile.jar
LEAD=<SET_VALUE> # value to trim from leading strand, often 20
TRAIN=<SET_VALUE> # value to trim from trailing strand, often 20
SLIDE=<SET_VALUE> # threshold and windlow length, often 4:20
MINREADLEN=<SET_VALUE> # minimum length for a read to be kept, often 90 for 150bp sequencing

# clipping
BAMUTILBAM=/path/to/bamutil/bin/bam/file

# bam statistics

# snp ID
ANGSD=~/programs/angsd/ # path to directory with angsd executables
SNPPVAL=<SET_VALUE> # max p-value for snp to be considered significant, often 1e-6
MINDEPTHIND=<SET_VALUE> # minimum depth per individual required for a site to be kept
MININD=<SET_VALUE> # minimum number of individuals required for a site to be kept
MINQ=<SET_VALUE> # minimum quality score required for a site to be kept
MINMAF=<SET_VALUE> # minimum minor allele frequency required for a site to be kept
MINMAPQ=<SET_VALUE> # minimum mapping quality score required for a site to be kept
POP=<SET_VALUE> # name of population