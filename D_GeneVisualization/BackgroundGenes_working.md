Given a list of sites with the first two columns describing 1-based chromosomal positions, filter out all sites corresponding to c!=3:
apply_mask_l mask_35_50.fa in.list > out.list  


grep 'ID=gene' genomic.gff > genomic.genes.gff

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





import pandas as pd

# Load the gene regions file (tab-separated)
gene_df = pd.read_csv("/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions.tsv", sep="\t") 
# Load the depth file (assumed to be tab-separated and without headers)
depth_df = pd.read_csv("/xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/all_avg_depthstats.txt", sep="\t", header=None, names=["Chromosome", "Position", "Depth"])

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
gene_df.to_csv("/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/gene_mask_proportions_depth.tsv", sep="\t", index=False)

#!/bin/bash

module load python
python makefile_gene_prop_depth.py


sbatch --account=mcnew \
        --job-name=submit_makefile_gene_prop_depth \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/submit_makefile_gene_prop_depth.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --mem=100gb \
        --time=10:00:00 \
        /xdisk/mcnew/dannyjackson/cardinals/datafiles/bamstats/submit_makefile_gene_prop_depth.sh