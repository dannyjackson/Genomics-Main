#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Command-line args: pre.mafs.gz post.mafs.gz
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript lrt_delta_af.R pre.mafs.gz post.mafs.gz")

dir <- args[1]
species <- args[2]

# Construct file paths
pre_file <- file.path(dir, paste0(species, "_pre.mafs.gz"))
post_file <- file.path(dir, paste0(species, "_post.mafs.gz"))

cat("Reading pre...\n")
pre <- fread(pre_file)

cat("Reading post...\n")
post <- fread(post_file)

cat("Merging...\n")
df <- merge(pre, post, by = c("chromo", "position"))

cat("Computing allele counts and LRT...\n")
df[, `:=`(
  n1 = 2 * nInd.x,
  n2 = 2 * nInd.y,
  k1 = round(2 * nInd.x * knownEM.x),
  k2 = round(2 * nInd.y * knownEM.y)
)]

# Shared freq under H0
df[, p_shared := (k1 + k2) / (n1 + n2)]

# Log likelihoods under H0 and H1
df[, logL0 := dbinom(k1, n1, p_shared, log = TRUE) + dbinom(k2, n2, p_shared, log = TRUE)]
df[, logL1 := dbinom(k1, n1, k1 / n1, log = TRUE) + dbinom(k2, n2, k2 / n2, log = TRUE)]

# LRT statistic and p-value
df[, LRT := 2 * (logL1 - logL0)]
df[, pval := pchisq(LRT, df = 1, lower.tail = FALSE)]

# Delta AF
df[, delta_af := knownEM.y - knownEM.x]

# Filter out sites where delta_af <= 0.25 before FDR correction
df_nonzero <- df[delta_af != 0 & abs(delta_af) > 0.25]

# FDR correction on nonzero delta AF sites
df_nonzero[, qval := p.adjust(pval, method = "fdr")]

# Save full output (with qval added only to filtered subset)
fwrite(df,  paste0(species, "/deltaAF_lrt_full.tsv"), sep = "\t")
fwrite(df_nonzero,  paste0(species, "/deltaAF_lrt_nonzero_deltaAF.tsv"), sep = "\t")

# Save filtered significant results
sig_p <- df_nonzero[pval < 0.001]
fwrite(sig_p,  paste0(species, "/deltaAF_lrt_significant_p0.001.tsv"), sep = "\t")

cat("Done. Significant SNPs (p < 0.001):", nrow(sig_p), "\n")

sig_p <- df_nonzero[pval < 0.01]
fwrite(sig_p,  paste0(species, "/deltaAF_lrt_significant_p0.01.tsv"), sep = "\t")

cat("Done. Significant SNPs (p < 0.01):", nrow(sig_p), "\n")

sig_p <- df_nonzero[pval < 0.05]
fwrite(sig_p,  paste0(species, "/deltaAF_lrt_significant_p0.05.tsv"), sep = "\t")

cat("Done. Significant SNPs (p < 0.05):", nrow(sig_p), "\n")

sig_q <- df_nonzero[qval < 0.20]
fwrite(sig_q,  paste0(species, "/deltaAF_lrt_significant_q0.20.tsv"), sep = "\t")

cat("Done. Significant SNPs (q < 0.20):", nrow(sig_q), "\n")

sig_q <- df_nonzero[qval < 0.10]
fwrite(sig_q,  paste0(species, "/deltaAF_lrt_significant_q0.10.tsv"), sep = "\t")

cat("Done. Significant SNPs (q < 0.10):", nrow(sig_q), "\n")

sig_q <- df_nonzero[qval < 0.05]
fwrite(sig_q,  paste0(species, "/deltaAF_lrt_significant_q0.05.tsv"), sep = "\t")

cat("Done. Significant SNPs (q < 0.05):", nrow(sig_q), "\n")