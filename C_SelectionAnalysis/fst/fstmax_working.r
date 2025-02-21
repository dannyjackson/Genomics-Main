#!/usr/bin/env Rscript

# max win fst 

# option 2
library(data.table)

# Define window size
window <- 25000  # Change this to your desired window size
win <- window/2  # Change this to your desired window size

# Read files efficiently
file1 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/25000/pyrrurban_pyrrrural.25000.fst", sep = "\t", header = TRUE)
file2 <- fread("/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/1/pyrrurban_pyrrrural.1.fst", sep = "\t", header = TRUE)


# Rename columns in file1 and file2 to avoid conflicts
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1"))
setnames(file2, c("region2", "chr2", "midPos2", "Nsites2", "fst2"))

# Set keys for fast searching
setkey(file2, chr2, midPos2)

# Function to find the maximum fst within the window
get_max_fst <- function(chr_val, pos_val) {
  max(file2[chr2 == chr_val & midPos2 >= (pos_val - win) & midPos2 <= (pos_val + win), fst2], na.rm = TRUE)
}

# Apply the function row-wise using a vectorized approach
file1[, fst_max := get_max_fst(chr1, midPos1), by = .(chr1, midPos1)]

# Rename columns back for output consistency
setnames(file1, c("region1", "chr1", "midPos1", "Nsites1", "fst1", "fst_max"),
                  c("region", "chr", "midPos", "Nsites", "fst", "fst_max"))

# Save the result
fwrite(file1, "/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.25000.fst.max", sep = "\t")

#!/bin/sh

module load R

Rscript /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/fstmaxwin.test.R



sbatch --account=mcnew \
        --job-name=fstmax.test \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/fstmax.test.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=3:00:00 \
        /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/fstmaxwin.test.sh




        # Local Score Approach
        # Load necessary libraries
library(dplyr)
library(ggplot2)

# Example: Simulated genomic scan data
set.seed(42)
genome_data <- data.frame(
  position = seq(1, 1e6, by = 1000),  # Positions every 1000 bp
  score = rnorm(1000, mean = 0, sd = 1) # Normally distributed random scores
)

# Introduce a selection signal (higher scores in a specific region)
genome_data$score[450:550] <- genome_data$score[450:550] + rnorm(101, mean = 3, sd = 0.5)

# Function to compute local score using cumulative sum
compute_local_score <- function(scores) {
  local_scores <- numeric(length(scores))
  max_score <- 0
  for (i in seq_along(scores)) {
    local_scores[i] <- max(0, local_scores[max(1, i - 1)] + scores[i])
    max_score <- max(max_score, local_scores[i])
  }
  return(local_scores)
}

# Compute local scores
genome_data$local_score <- compute_local_score(genome_data$score)

# Permutation test to determine significance threshold
num_permutations <- 1000
max_null_scores <- numeric(num_permutations)

for (perm in 1:num_permutations) {
  shuffled_scores <- sample(genome_data$score)  # Shuffle scores
  max_null_scores[perm] <- max(compute_local_score(shuffled_scores))
}

# Determine empirical significance threshold (e.g., 95th percentile)
threshold <- quantile(max_null_scores, 0.95)

# Identify significant regions
genome_data$significant <- genome_data$local_score > threshold

# Plot results
ggplot(genome_data, aes(x = position, y = local_score)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_point(data = genome_data[genome_data$significant, ], aes(x = position, y = local_score), color = "red") +
  theme_minimal() +
  labs(title = "Local Score Approach in Genome Scans",
       x = "Genomic Position",
       y = "Local Score",
       caption = "Dashed line: 95% significance threshold")
