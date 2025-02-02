#!/usr/bin/env Rscript

# Load required packages, installing if necessary
required_packages <- c("qqman", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer")
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
pop1 <- args[2]
pop2 <- args[3]
color1 <- args[4]
color2 <- args[5]
cutoff <- as.numeric(args[6])  # Convert to numeric
input <- args[7]

# Detect file type based on header
cat("Detecting input data type...\n")
data <- read.delim(input, sep = '\t') %>% na.omit()


if ("dxy" %in% names(data)) {
  metric <- "dxy"
} else if ("fst" %in% names(data)) {
  metric <- "fst"
} else if ("tajimas_d" %in% names(data)) {
  metric <- "tajimas_d"
} else {
  stop("Unknown data format. Ensure the input file contains dxy, fst, or tajima's D column.")
}

# Z-transform metric values
cat("Calculating Z-transform and identifying top outliers...\n")
metric_xbar <- mean(data[[metric]], na.rm = TRUE)
metric_sd <- sd(data[[metric]], na.rm = TRUE)
data$z <- (data[[metric]] - metric_xbar) / metric_sd
data$neg_log_pvalues_one_tailed <- -log10(pnorm(data$z, lower.tail = FALSE))

# Identify top outliers
ntotal <- nrow(data)
top_snps <- round(ntotal * cutoff)
outlier_data <- data %>% 
  arrange(desc(neg_log_pvalues_one_tailed)) %>% 
  head(top_snps) %>% 
  arrange(chromo, position)

metric_cutoff <- min(outlier_data[[metric]])

# Save cutoff value
cat("Saving cutoff value...\n")
cutoff_file <- file.path(outdir, "analyses", metric, paste0(pop1, "_", pop2, "_", metric, "_stats.txt"))
cat(metric, "cutoff:", metric_cutoff, "\n", file = cutoff_file, append = TRUE)

# Save outliers
cat("Saving outliers data...\n")
outlier_file <- file.path(outdir, "analyses", metric, paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".", metric, ".snps.outlier.csv"))
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
  geom_hex(bins = 100) +  # Fast binning for dense regions
  geom_point(data = outlier_dxy, aes(x = BPcum, y = dxy)), +
  scale_fill_viridis_c() +  # Heatmap coloring for density
  scale_color_manual(values = rep(c(color1, color2), length(unique(data$chromo)) / 2)) +
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

ggsave(filename = file.path(outdir, "analyses", metric, paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".", metric, ".snps.sigline.png")), 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
