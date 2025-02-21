# SNPeff

# install newer java version
wget https://download.oracle.com/java/23/latest/jdk-23_linux-x64_bin.tar.gz
gunzip jdk-23_linux-x64_bin.tar.gz

# module load openjdk/19.0.1 ... must be 21 or newer

### snpeff
# https://hackmd.io/@tlama/outlierFST#Evaluate-outliers-in-R-see-Rmarkdown-here
# downloaded snpeff
# download snpeff databases for chicken and great tit
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

git clone https://github.com/pcingola/SnpEff.git
module load maven/3.6.3

mvn install:install-file \
    -Dfile=biojava3-core-3.0.7.jar \
    -DgroupId=org.biojava \
    -DartifactId=biojava3-core \
    -Dversion=3.0.7 \
    -Dpackaging=jar

# Annotate SNPs in SnpEff (Cingolani et al. 2012)
# check if prebuilt genome exists in database
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar databases | grep -i 'Camarhynchus'
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar databases | grep -i 'parvulus'
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar databases | grep -i 'finch'
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar databases | grep -i 'GCF_901933205'

# nope so we have to build it
# added these two lines to snpEff.config:

# Small tree finch genome, version GCA_901933205
GCA_901933205.genome : SmallTreeFinch

# Get the genome and uncompress it:

# Create directory for this new genome
cd /home/u15/dannyjackson/programs/snpEff/data
mkdir GCA_901933205
cd GCA_901933205
# Get annotation files
cp /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff ./genes.gff
mkdir genomes
cd genomes
cp /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna ./GCA_901933205.fa
cp /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/protein.faa ./protein.fa
cp /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/cds_from_genomic.fna ./cds.fa


sed -i '$ d' snpEff.config

echo 'GCA_901933205.genome : GCA_901933205' >> snpEff.config

~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff3 -v GCA_901933205 2>&1 | tee GCA_901933205.build.output

sed -i 's/java/\~\/programs\/jdk\-23\.0\.2\/bin\/java/g' /home/u15/dannyjackson/programs/snpEff/scripts/buildDbNcbi.sh 

sed -i 's/jdk\-2/jdk\-23/g' /home/u15/dannyjackson/programs/snpEff/scripts/buildDbNcbi.sh 

sed -i '$ d' snpEff.config

/home/u15/dannyjackson/programs/snpEff/scripts/buildDbNcbi.sh GCF_901933205.1

FATAL ERROR: No CDS checked. This might be caused by differences in FASTA file transcript IDs respect to database's transcript's IDs.
Transcript IDs from database (sample):
        'XM_030965427.1'
        'rna-XM_030956335.1'
        'XM_030965425.1'
        'XM_030965424.1'

Transcript IDs from database (fasta file):
        'lcl|NC_044572.1_cds_XP_030823689.1_2603'
        'lcl|NC_044586.1_cds_XP_030815405.1_19420'
        'lcl|NC_044571.1_cds_XP_030801370.1_1690'
        'lcl|NC_044575.1_cds_XP_030805537.1_8939'


# fix cds
import sys

def restructure_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.split()
                protein_id = None
                for part in parts:
                    if part.startswith("[protein_id="):
                        protein_id = part.split("=")[1].rstrip("]")
                        break
                if protein_id:
                    new_header = f">{protein_id} {line[1:]}"
                    outfile.write(new_header)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

if __name__ == "__main__":
    input_fasta = "cds.fa"  # Change to your input file name
    output_fasta = "cds.revised.fa"  # Change to your desired output file name
    restructure_fasta(input_fasta, output_fasta)



# this improved the situation but I am now getting this error:
FATAL ERROR: No CDS checked. This might be caused by differences in FASTA file transcript IDs respect to database's transcript's IDs.
Transcript IDs from database (sample):
        'XM_030965427.1'
        'rna-XM_030956335.1'
        'XM_030965425.1'
        'XM_030965424.1'
        'rna-XR_004060613.1'
        'rna-XM_030952858.1'
Transcript IDs from database (fasta file):
        'XP_030820196.1'
        'XP_030822353.1'
        'XP_030800543.1'


# fix gtf
import re

def replace_transcript_with_protein(gtf_file, output_file):
    # Open the input GTF file and the output file
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip comment lines
            if line.startswith("#") or line.startswith("!"):
                outfile.write(line)
                continue

            # Match and extract transcript_id and protein_id
            transcript_match = re.search(r'transcript_id\s+"([^"]+)"', line)
            protein_match = re.search(r'protein_id\s+"([^"]+)"', line)

            # If both transcript_id and protein_id are found, replace transcript_id with protein_id
            if transcript_match and protein_match:
                transcript_id = transcript_match.group(1)
                protein_id = protein_match.group(1)
                # Replace transcript_id with protein_id in the line
                line = line.replace(f'transcript_id "{transcript_id}"', f'transcript_id "{protein_id}"')

            # Write the modified (or unmodified) line to the output file
            outfile.write(line)

# Example usage
input_gtf_file = '../genes.gtf'
output_gtf_file = 'genes.gtf'
replace_transcript_with_protein(input_gtf_file, output_gtf_file)



~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff3 -v GCA_901933205 2>&1 | tee GCA_901933205.build.output
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff2 -v GCA_901933205 2>&1 | tee GCA_901933205.build.output

cd /home/u15/dannyjackson/programs/snpEff
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v GCA_901933205 2>&1 | tee GCA_901933205.build.output
# broke on memory

interactive -a mcnew -m 50gb
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v GCA_901933205 2>&1 | tee GCA_901933205.build.output





WARNING_TRANSCRIPT_NOT_FOUND: Cannot find transcript 'XM_030955085.1'. Created transcript 'XM_030955085.1' and gene 'LOC115907285' for NC_044571.1      Gnomon  EXON    1693    3360
    -
        db_xref : GeneID:115907285
        exon_number : 1
        gbkey : mRNA
        gene : LOC115907285
        gene_id : LOC115907285
        model_evidence : Supporting evidence includes similarity to: 95% coverage of the annotated genomic feature by RNAseq alignments, including 3 samples with support for all annotated introns
        product : uncharacterized LOC115907285
        transcript_id : XM_030955085.1
. File '/home/u15/dannyjackson/programs/snpEff/./data/GCA_901933205/genes.gtf' line 6   'NC_044571.1    Gnomon  exon    1694    3361    .       -       .       gene_id "LOC115907285"; transcript_id "XM_030955085.1"; db_xref "GeneID:115907285"; gbkey "mRNA"; gene "LOC115907285"; model_evidence "Supporting evidence includes similarity to: 95% coverage of the annotated genomic feature by RNAseq alignments, including 3 samples with support for all annotated introns"; product "uncharacterized LOC115907285"; exon_number "1"; '


# jfc that was awful but now it works
# i don't know what's up with these transcript issues but I'm choosing to ignore it

# annotate all snps in vcf

#!/bin/sh

~/programs/jdk-23.0.2/bin/java -jar /home/u15/dannyjackson/programs/snpEff/snpEff.jar -c /home/u15/dannyjackson/programs/snpEff/snpEff.config -v GCA_901933205 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.vcf > /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.snpEffann.vcf 

grep 'missense_variant' /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.snpEffann.vcf 

sbatch --account=mcnew \
        --job-name=snpeff \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.snpeff.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=1:00:00 \
        --mem=50gb \
        /xdisk/mcnew/dannyjackson/cardinals/scripts/snpeff_wholegenome.sh

#!/bin/sh

~/programs/jdk-23.0.2/bin/java -jar /home/u15/dannyjackson/programs/snpEff/snpEff.jar -c /home/u15/dannyjackson/programs/snpEff/snpEff.config -v GCA_901933205 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf 

sbatch --account=mcnew \
        --job-name=snpeff \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.snpeff.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=5:00:00 \
        --mem=50gb \
        /xdisk/mcnew/dannyjackson/cardinals/scripts/snpeff_filtered.sh

# now filter fst file to contain only genes within snpeff vcf

