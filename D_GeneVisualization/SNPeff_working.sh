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

cd /home/u15/dannyjackson/programs/snpEff/

sed -i '$ d' snpEff.config

echo 'GCA_901933205.genome : GCA_901933205' >> snpEff.config
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff3 -v GCA_901933205
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -genbank GCF_901933205

# annotate all snps in vcf
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar -c /home/u15/dannyjackson/programs/snpEff/snpEff.config -v GCA_901933205 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.vcf > /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.snpEffann.vcf 

~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff3 -v ncbi_cardinalis 
~/programs/jdk-23.0.2/bin/java -Xmx4g -jar /path/to/snpEff/snpEff.jar -c /path/to/snpEff/snpEff.config -v b10k VYXE01005221.vcf > VYXE01005221.an…
https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_b10k.zip











# I need to make a vcf of outlier snps
# in this file, $2 is chromo and $3 is position
# /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.Ztransformed.csv
# I need to convert chromosome names in the output file to new chrom names 

CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
FILE="/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.Ztransformed.csv"
OUT_FILE="/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.Ztransformed.chrnames.csv"

# Process the CSV file, replacing only the second column
awk -F',' -v CHR_FILE="$CHR_FILE" '
BEGIN {
    OFS = "\t";  # Set output field separator to tab
    while ((getline < CHR_FILE) > 0) {
        map[$1] = $2;  # Store mappings from second column to first column
    }
}
BEGIN { FS = "\t" }  # Set input field separator to tab for FILE processing
{
    if ($2 in map) {
        $2 = map[$2];  # Replace second column if match found
    }
    print;
}' "$FILE" > "$OUT_FILE"

echo "Finished processing. Output saved to $OUT_FILE"

# make a regions file from this
# CHROM, BEG, END
# given a file with chromo in $2 and position in $3,
awk '{print $2, $3}' /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.Ztransformed.chrnames.csv | awk 'BEGIN {OFS="\t"} NR==1 {print "CHROM", "BEG", "END"; next} {print $1, $2-12500, $2+12500}' > /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.outlier.regionsfile.tsv

bcftools view [OPTIONS] file.vcf.gz [REGION […​]]


bcftools view /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_snps_multiallelic.vcf -R /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/nocaurban_nocarural.fst.25000.outlier.regionsfile.tsv -o nocaurban_nocarural.fst.25000.outliers.vcf

# ~/programs/jdk-23.0.2/bin/java -jar snpEff.jar download GRCg6a.105
~/programs/jdk-23.0.2/bin/java -Xmx20g -jar snpEff.jar build -v b10k 

~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v b10k -noCheckProtein
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v b10k
~/programs/jdk-23.0.2/bin/java -Xmx8g -jar snpEff.jar -v b10k VYXE01005221.vcf > VYXE01005221.ann.vcf
~/programs/jdk-23.0.2/bin/java -Xmx8g -jar snpEff.jar -v b10k noca.filtered.geno25.maf1.vcf > noca.filtered.geno25.maf1.ann.vcf
~/programs/jdk-23.0.2/bin/java -Xmx8g -jar snpEff.jar -v b10k pyrr.filtered.geno25.maf1.vcf > pyrr.filtered.geno25.maf1.ann.vcf
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gff3 -v ncbi_cardinalis 
~/programs/jdk-23.0.2/bin/java -Xmx4g -jar /path/to/snpEff/snpEff.jar -c /path/to/snpEff/snpEff.config -v b10k VYXE01005221.vcf > VYXE01005221.an…
https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_b10k.zip
. File '/Users/danjack/snpEff/./data/b10k/genes.gtf' line 33	'VYXE01000004.1	Genbank	CDS	19854	19961	.	-	0	gene_id	"CAR…
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v b10k
~/programs/jdk-23.0.2/bin/java -jar snpEff.jar build -gtf22 -v b10k -noCheckProtein