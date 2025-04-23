#!/usr/bin/env Rscript

# Load required library
suppressPackageStartupMessages(library(data.table))

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compute_delta_af.R <dir> <species>")
}

dir <- args[1]
species <- args[2]

# Construct file paths
pre_file <- file.path(dir, paste0(species, "_pre.mafs"))
post_file <- file.path(dir, paste0(species, "_post.mafs"))

# Read only needed columns
cat("Reading pre file...\n")
pre <- fread(pre_file, select = c("chromo", "position", "major", "minor", "knownEM", "nInd"))

cat("Reading post file...\n")
post <- fread(post_file, select = c("chromo", "position", "major", "minor", "knownEM", "nInd"))

# Merge by chromosome and position
cat("Merging files...\n")
merged <- merge(pre, post, by = c("chromo", "position"))

# Check for mismatches in major and minor alleles across time points
cat("Checking allele identity consistency...\n")
major_mismatch <- merged[major.x != major.y]
minor_mismatch <- merged[minor.x != minor.y]

cat("Major allele mismatches:", nrow(major_mismatch), "\n")
cat("Minor allele mismatches:", nrow(minor_mismatch), "\n")

if (nrow(major_mismatch) > 0 || nrow(minor_mismatch) > 0) {
  cat("Warning: some sites have flipped major/minor alleles.\n")
  cat("Consider resolving these before interpreting delta AF.\n")
} else {
  cat("All major/minor alleles match across timepoints. ✅\n")
}

# Compute ΔAF
cat("Computing delta allele frequency...\n")
merged[, delta_af := knownEM.y - knownEM.x]
merged <- merged[!is.na(delta_af)]

# Estimate SE using binomial approximation (assuming HWE and diploid)
merged[, pre_se := sqrt(knownEM.x * (1 - knownEM.x) / (2 * nInd.x))]
merged[, post_se := sqrt(knownEM.y * (1 - knownEM.y) / (2 * nInd.y))]
merged[, delta_se := sqrt(pre_se^2 + post_se^2)]

# Z-test
merged[, z := delta_af / delta_se]
merged[, p := 2 * pnorm(-abs(z))]
merged[, q := p.adjust(p, method = "fdr")]

# Filter for significant sites
significant_sites <- merged[q < 0.05 & !is.na(q)]

# Output
fwrite(merged, file.path(dir, paste0(species, "_deltaAF_with_significance.tsv")), sep = "\t")
fwrite(significant_sites, file.path(dir, paste0(species, "_sig_deltaAF.tsv")), sep = "\t")

cat("Done.\n")
