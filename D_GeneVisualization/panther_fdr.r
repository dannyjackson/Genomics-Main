#!/usr/bin/env Rscript

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

df <- read.csv(paste0(file, ".txt"), sep = '\t')
pvals <- df$P.value

qvals <- p.adjust(pvals, method = "fdr")

df$qvals <- qvals

df_filtered = df[qvals < 0.2 ,]

write.csv(df_filtered, paste0(file, ".fdr.txt"), quote = FALSE, row.names = FALSE)
