#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
stats_file <- paste0(output, ".stats")

# Load data
df <- read.table(input, header=FALSE, sep="\t")

# Assign column names
colnames(df) <- c("chrom", "win_start", "win_end", "avg_stat", "numsites")

# Calculate mean and standard deviation
mean_depth <- mean(df$avg_stat, na.rm=TRUE)
sd_depth <- sd(df$avg_stat, na.rm=TRUE)

# Define cutoffs (±2 SD)
lower_cutoff <- mean_depth - 2 * sd_depthcd
upper_cutoff <- mean_depth + 2 * sd_depth

# Filter data
filtered_df <- df[df$avg_stat >= lower_cutoff & df$avg_stat <= upper_cutoff, c(1,2,3)]

# Record changes in file size
df_length <- nrow(df)
filtered_df_length <- nrow(filtered_df)
diff_length <- df_length - filtered_df_length

# Save filtered data
write.table(filtered_df, output, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

# Print cutoffs
writeLines(c(
  paste("Lower cutoff:", lower_cutoff),
  paste("Upper cutoff:", upper_cutoff),
  paste("Original file length:", df_length),
  paste("Filtered file length:", filtered_df_length),
  paste("Number of rows removes:", diff_length)
), con = stats_file)
