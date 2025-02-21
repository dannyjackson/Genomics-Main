Given a list of sites with the first two columns describing 1-based chromosomal positions, filter out all sites corresponding to c!=3:
apply_mask_l mask_35_50.fa in.list > out.list  

# start by just generating a basic list of all genes in the reference genome
grep 'ID=gene' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff > /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff

awk '{print $9}' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff | awk -F"[-;]" '{print $2}' | sort -u > genelist.txt

import pysam
import csv

# Input file paths
gene_file = "/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff"  # Assumed tab-separated
mask_fasta = "/xdisk/mcnew/dannyjackson/cardinals/datafiles/snpable/GCF_901933205_mask.150.50.fa"

# Input file paths

output_file = "/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv"

# Load the mask FASTA file
mask = pysam.FastaFile(mask_fasta)

# Process the gene file and write to output
with open(gene_file, "r") as f, open(output_file, "w", newline="") as out_f:
    reader = csv.reader(f, delimiter="\t")
    writer = csv.writer(out_f, delimiter="\t")
    # Write header
    writer.writerow(["Chromosome", "Start", "End", "Proportion_3", "Gene"])
    for row in reader:
        chrom, start, end, gene = row[0], int(row[3]), int(row[4]), row[8]  # Extract chromosome, start, end
        # Extract the masked sequence for this region
        sequence = mask.fetch(chrom, start - 1, end)  # Convert to 0-based indexing
        # Count occurrences of "3"
        count_3 = sequence.count("3")
        total_length = len(sequence)
        proportion = count_3 / total_length if total_length > 0 else 0
        # Write results to file
        writer.writerow([chrom, start, end, proportion, gene])

print(f"Results saved to {output_file}")


awk -F"\t" '$4 > 0.99' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv | wc -l 
15151/16831 = 90% of genes > 0.99

awk -F"\t" '$4 > 0.95' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv | wc -l 
15921/16831 = 94.6% of genes > 0.95


# create column with average depth per gene
# first, create a file of depth for all individuals
#!/bin/bash
module load samtools

samtools depth /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/*.realigned.bam >> /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_depthstats.txt 

sbatch --account=mcnew \
        --job-name=depthstats \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/depthstats.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=3:00:00 \
        /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/computeavgdepthstats.sh

# compute average depth
#!/bin/bash

awk 'NR==1 {next} { 
    sum = 0; 
    for (i = 3; i <= NF; i++) sum += $i; 
    avg = sum / (NF - 2); 
    print $1, $2, avg;
}' /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_depthstats.txt  > /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.txt 

sbatch --account=mcnew \
        --job-name=depthstats \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/depthstats.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=3:00:00 \
        /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/avgdepthstats.sh





import sys
import pandas as pd

s = sys.argv[1]  # Get the first command-line argument
gene_file_path = f"/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.{s}.tsv"

gene_df = pd.read_csv(gene_file_path, sep="\t")

depth_file_path = f"/xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.{s}.tsv"

depth_df = pd.read_csv(depth_file_path, sep="\t", header=None, names=["Chromosome", "Position", "Depth"])

# Function to compute average depth for each gene
def compute_avg_depth(row):
    subset = depth_df[
        (depth_df["Chromosome"] == row["Chromosome"]) &
        (depth_df["Position"] >= row["Start"]) &
        (depth_df["Position"] <= row["End"])
    ]
    return subset["Depth"].mean() if not subset.empty else 0  # Return 0 if no depths found

# Apply the function to each row in the gene DataFrame
gene_df["Average_Depth"] = gene_df.apply(compute_avg_depth, axis=1)
# Save the result
output_file_path = f"/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.{s}.tsv"

gene_df.to_csv(output_file_path, sep="\t", index=False)


for s in `cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt`;
	do echo $s

    grep $s /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv > "/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.$s.tsv"

    grep $s /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.txt > "/xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.$s.txt"

done

tail -1 /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt > /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.test.txt

for s in `cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt`;
	do echo ${s}

	sbatch --account=mcnew \
	--job-name=submit_makefile_gene_prop_depth_${s} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.submit_makefile_gene_prop_depth_${s}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
    --mem=50gb \
	--time=6:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/submit_makefile_gene_prop_depth.sh -s "$s"

done

# 12228322

#!/bin/bash

# Parse command-line arguments
while getopts "s:" option; do
    case "${option}" in
        s) SCAFF=${OPTARG} ;;
        *) usage ;;
    esac
done

module load python

echo ${SCAFF}

echo "first file"
head -1 /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv > "/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.${SCAFF}.tsv"
grep ${SCAFF} /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv >> "/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.${SCAFF}.tsv"

echo "second file"
grep ${SCAFF} /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.txt > "/xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.${SCAFF}.txt"

echo "python"
python /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/makefile_gene_prop_depth.py ${SCAFF}



cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/

head -n 1 gene_mask_proportions_depth.NC_044571.1.tsv > gene_mask_proportions_depth.all.txt

for s in `cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/SCAFFOLDS.txt`;
	do echo ${s}

	cat "gene_mask_proportions_depth.$s.tsv" | tail -n +2 >> gene_mask_proportions_depth.all.txt

done

# Make histograms
cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats

# first, identify genes with high depth:
28.15031

head -n 1 /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.txt > /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.filtereddepthmax.txt

awk '$6 < 28.15031 { print $0 }' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.txt >> /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.filtereddepthmax.txt

# cutoff of 2xSD is 28.15031
# cutoff of 3xSD is 39.76541
# a bunch of LOC genes and also:
ACAD9       50.778   Acyl-CoA dehydrogenase
# ARF4        32.566   Auxin response factors (ARF) GTPase
CEP250      69.202 Centrosomal Protein 250
# GARNL3      35.596
# MMP17       32.290
# NAA20       33.701
PROCA1      70.155
RESF1       43.456
# TRNAQ-CUG   32.025
# UBE2G2      30.689

# go with 2x...


# what genes had low mapability?
awk '$4 < 0.5 { print $0 }' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.txt | wc -l

awk '($4 < 0.95) { match($5, /ID=gene-([^;]+)/, arr); print arr[1] }' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.txt | grep -v 'LOC'

R

library(ggplot2)

# Read your data
df <- read.csv("/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.filtereddepthmax.txt", sep='\t')

df$Average_Depth <- as.numeric(df$Average_Depth)

df_plot = df[!is.na(df$Average_Depth),]


u = mean(df_plot$Average_Depth)
# 4.568
stdv = sd(df_plot$Average_Depth)
# 1.221131
u + (2*stdv)
# cutoff of 28.15031 depth originally
# 7.010419 after filtering



# Histogram for Depth
p2 <- ggplot(df, aes(x = Average_Depth)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  theme_minimal() +
  ggtitle("Histogram of Depth")

# Save the Depth histogram
ggsave("Depth_Histogram.png", plot = p2, width = 6, height = 4, dpi = 300)






df$Proportion_3 <- as.numeric(df$Proportion_3)

# Histogram for Proportion
p1 <- ggplot(df, aes(x = Proportion_3)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  theme_minimal() +
  ggtitle("Histogram of Proportion")

# Save the Proportion histogram
ggsave("Proportion_Histogram.png", plot = p1, width = 6, height = 4, dpi = 300)


# decide to filter based on 3 < depth < 28.15031 and proportion > 0.95, see what genes that cuts out

awk '($6 < 3 || $6 > 28.15031 || $4 < 0.90)  { match($5, /ID=gene-([^;]+)/, arr); print arr[1] }'  /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.all.txt | grep -v 'LOC' | sort -u > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/excludedgenes.3depth28.prop90.txt

