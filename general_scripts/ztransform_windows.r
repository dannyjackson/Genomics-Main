#!/usr/bin/env Rscript

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table", "tidyr")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
cutoff <- as.numeric(args[2])  # Convert to numeric
input <- args[3]
win <- args[4]
pop1 <- args[5]
pop2 <- ifelse(length(args) > 5 && args[6] != "", args[6], NA)

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Detect file type based on header
cat("Detecting input data type...\n")

data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

new_names <- names(data) %>%
  gsub("(?i)\\bU\\b", "raisd", ., perl = TRUE)

names(data) <- new_names

cat("read in input data...\n")

if ("dxy" %in% names(data)) {
  metric <- "dxy"
} else if ("fst" %in% names(data)) {
  metric <- "fst"
} else if ("fst_max" %in% names(data)) {
  metric <- "fst_max"
} else if ("Tajima" %in% names(data)) {
  metric <- "Tajima"
} else if ("raisd" %in% names(data)) {
  metric <- "raisd"
} else if ("delta_af" %in% names(data)) {
  metric <- "delta_af"
} else {
  stop("Unknown data format. Ensure the input file contains dxy, fst, or tajima's D column.")
}

# standardize use of chromosome and pos in file

cat("renaming names...\n")

# Define patterns and replacements
new_names <- names(data) %>%
  gsub("(?i)\\bchr(?:om(?:osome)?)?\\b", "chromo", ., perl = TRUE) %>%   # Replace chr, chrom, and chromosome with "chromo"
  gsub("(?i)\\b(?:mid|pos|midpos|WinCenter)\\b", "position", ., perl = TRUE)  # Replace pos, mid, midpos variants with "position"


# Assign new column names
names(data) <- new_names

# option 1
# Z-transform metric values

# Compute mean and SD 
cat("Computing mean and SD...\n")
metric_xbar <- mean(data[[metric]], na.rm = TRUE)
metric_sd <- sd(data[[metric]], na.rm = TRUE)

cat("Calculating Z-transform and identifying top outliers...\n")

metric_xbar <- mean(data[[metric]], na.rm = TRUE)
metric_sd <- sd(data[[metric]], na.rm = TRUE)
data$z <- (data[[metric]] - metric_xbar) / metric_sd
data$neg_log_pvalues_one_tailed <- -log10(pnorm(data$z, lower.tail = FALSE))

df <- data[ -c(1) ]

# save file
cat("Saving Z-transformed data...\n")
z_file <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", pop_name, ".", metric, ".", win, ".Ztransformed.csv"))
write_delim(df, z_file, quote = "none", delim = "\t")
  