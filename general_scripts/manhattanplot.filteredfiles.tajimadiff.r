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
metric <- args[6]

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data[!is.na(Tajima)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(Tajima)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff_top <- min(top_snps_dt[[metric]], na.rm = TRUE)


# Save cutoff value
cat("Saving cutoff value...\n")
cutoff_file <- paste0(input, ".top.stats.txt")
cat(metric, "cutoff:", metric_cutoff_top, "\n", file = cutoff_file, append = TRUE)

# Save outliers
cat("Saving outliers data...\n")
outlier_file <- paste0(input, ".top.outlier.csv")
write.csv(top_snps_dt, outlier_file, row.names = FALSE)

# identify bottom snps

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data[!is.na(Tajima)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(Tajima)) %>%
  slice_tail(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff_bottom <- min(top_snps_dt[[metric]], na.rm = TRUE)


# Save cutoff value
cat("Saving cutoff value...\n")
cutoff_file <- paste0(input, ".bottom.stats.txt")
cat(metric, "cutoff:", metric_cutoff_bottom, "\n", file = cutoff_file, append = TRUE)

# Save outliers
cat("Saving outliers data...\n")
outlier_file <- paste0(input, ".bottom.outlier.csv")
write.csv(top_snps_dt, outlier_file, row.names = FALSE)

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
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff_top) +
  geom_hline(yintercept = metric_cutoff_bottom) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = "${input}.sigline.png", 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
