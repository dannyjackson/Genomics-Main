#!/usr/bin/env Rscript

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
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
color1 <- args[2]
color2 <- args[3]
cutoff <- as.numeric(args[4])  # Convert to numeric
input <- args[5]
win <- args[6]
pop1 <- args[7]
pop2 <- ifelse(length(args) > 7 && args[8] != "", args[8], NA)

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Detect file type based on header
cat("Detecting input data type...\n")

data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("read in input data...\n")

if ("dxy" %in% names(data)) {
  metric <- "dxy"
} else if ("fst" %in% names(data)) {
  metric <- "fst"
} else if ("Tajima" %in% names(data)) {
  metric <- "Tajima"
} else {
  stop("Unknown data format. Ensure the input file contains dxy, fst, or tajima's D column.")
}

# standardize use of chromosome and pos in file

cat("renaming names...\n")

# Define patterns and replacements
new_names <- names(data) %>%
  gsub("(?i)\\bchr(?:omosome)?\\b", "chromo", ., perl = TRUE) %>%   # Replace chromosome variants with "chromo"
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

# Compute z-scores only when needed (avoiding extra column storage)
chunk_size <- 1e6  # Adjust based on available memory
num_rows <- nrow(data)
num_chunks <- ceiling(num_rows / chunk_size)

for (i in seq_len(num_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, num_rows)
  
  data[start:end, neg_log_pvalues_one_tailed := 
    -log10(pnorm((get(metric) - metric_xbar) / metric_sd, lower.tail = FALSE))]
  
  # Print an update every 100 chunks
  if (i %% 100 == 0 || i == num_chunks) {
    cat(sprintf("Processed %d/%d chunks (%.2f%% complete)\n", i, num_chunks, (i / num_chunks) * 100))
    flush.console()  # Ensure the output is printed immediately
  }
}

# Identify top outliers
ntotal <- nrow(data)
top_snps <- round(ntotal * cutoff)

outlier_data <- data[order(-neg_log_pvalues_one_tailed)][1:top_snps][order(chromo, position)]

metric_cutoff <- min(outlier_data[[metric]], na.rm = TRUE)







# Save cutoff value
cat("Saving cutoff value...\n")
cutoff_file <- file.path(outdir, "analyses", metric, paste0(pop_name, "_", metric, "_", win, "_stats.txt"))
cat(metric, "cutoff:", metric_cutoff, "\n", file = cutoff_file, append = TRUE)

# Save outliers
cat("Saving outliers data...\n")
outlier_file <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", pop_name, ".", metric, "_", win, ".snps.outlier.csv"))
write.csv(outlier_data, outlier_file, row.names = FALSE)

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chromo <- factor(data$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")
ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = file.path(outdir, "analyses", metric, paste0(pop_name, "/", win, "/", pop_name, ".", metric, ".", win, ".sigline.png")), 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
