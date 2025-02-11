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
metric <- args[7]
pop1 <- args[8]
pop2 <- ifelse(length(args) > 8 && args[9] != "", args[9], NA)

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data <- fread(input, sep = ",", data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
top_snps_count <- round(nrow(data) * cutoff)
cat("identifying top snps 2...\n")
top_snps_dt <- data[order(-neg_log_pvalues_one_tailed, na.last = TRUE)][1:top_snps_count]

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- top_snps_dt[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)





# Save cutoff value
cat("Saving cutoff value...\n")
cutoff_file <- file.path(outdir, "analyses", metric, paste0(pop_name, "_", metric, "_", win, "_stats.txt"))
cat(metric, "cutoff:", metric_cutoff, "\n", file = cutoff_file, append = TRUE)

# Save outliers
cat("Saving outliers data...\n")
outlier_file <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", pop_name, ".", metric, "_", win, ".outlier.csv"))
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
