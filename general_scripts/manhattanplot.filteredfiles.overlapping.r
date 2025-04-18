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

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
color1 <- args[1]
color2 <- args[1]
cutoff <- as.numeric(args[1])  # Convert to numeric
input1 <- args[1]
input2 <- args[1]
metric <- "raisd"
outfile <- args[1]

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)



# Prepare data for plotting
# data 1
cat("Preparing data for plotting...\n")
data1$chromo <- factor(data1$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data1 <- data1 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data1, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data1 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# data 2
cat("Preparing data for plotting...\n")
data2$chromo <- factor(data2$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data2 <- data2 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf2 <- plot_data2 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)

  geom_point(data = plot_data1, aes(x = BPcum, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = plot_data2, aes(x = BPcum, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data1$chromo)))) +
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


ggsave(filename = outfile, 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
