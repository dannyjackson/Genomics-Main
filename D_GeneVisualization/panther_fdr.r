#!/usr/bin/env Rscript

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

df <- read.csv(paste0(file, ".txt"), sep = '\t')
pvals <- df$P.value

qvals <- p.adjust(pvals, method = "fdr")

df$qvals <- qvals

df_filtered_05 = df[qvals < 0.05 ,]

write.csv(df_filtered_05, paste0(file, ".fdr.05.txt"), quote = FALSE, row.names = FALSE)

df_filtered_1 = df[qvals < 0.1 ,]

write.csv(df_filtered_1, paste0(file, ".fdr.1.txt"), quote = FALSE, row.names = FALSE)

df_filtered_2 = df[qvals < 0.2 ,]

write.csv(df_filtered_2, paste0(file, ".fdr.2.txt"), quote = FALSE, row.names = FALSE)
