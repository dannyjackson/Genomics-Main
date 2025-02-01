#!/usr/bin/env Rscript
# Load required packages, installing if necessary
required_packages <- c("qqman", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer")
installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pop1 <- args[2]
pop2 <- args[3]
color1 <- args[4]
color2 <- args[5]
cutoff <- args[6]

# Read and clean data
dxy_file <- file.path(outdir, "analyses/dxy", paste0(pop1, "_", pop2, "/Dxy_persite_", pop1, "_", pop2, ".autosomes.txt"))
dxy <- read.csv(dxy_file, sep = '\t') %>% na.omit()

dxy_file <-
  read.csv(paste0(outdir, "/analyses/dxy/",
                  pop1, "_", pop2,
                  "/Dxy_persite_", pop1, "_", pop2, ".autosomes.txt"),
                  sep = "\t")  %>% na.omit()

min_dxy <- min(dxy$dxy)
max_dxy <- max(dxy$dxy)

cat(c("Min dxy:", min_dxy),
    file = paste0(outdir, "/analyses/dxy/",pop1, "_", pop2, "/", pop1, "_", pop2, "dxy_stats.txt"),
    sep = "\n", append = TRUE)

cat(c("Max dxy:", max_dxy),
    file = paste0(outdir, "/analyses/dxy/",pop1, "_", pop2, "/", pop1, "_", pop2, "dxy_stats.txt"),
    sep = "\n", append = TRUE)


# Z-transform dxy values
dxy_xbar <- mean(dxy$dxy, na.rm = TRUE)
dxy_sd <- sd(dxy$dxy, na.rm = TRUE)
dxy$z <- (dxy$dxy - dxy_xbar) / dxy_sd
p_values_one_tailed <- pnorm(q = dxy$z, lower.tail = FALSE)

# Calculate -log10 of the p-value
dxy$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

# Identify top outliers
nsnps <- nrow(dxy)
top_snps <- round(nsnps * cutoff)
outlier_dxy <- dxy %>% 
  arrange(desc(neg_log_pvalues_one_tailed)) %>% 
  head(top_snps) %>% 
  arrange(chromo, position)

dxy_cutoff <- min(outlier_dxy$dxy)

# Save cutoff value
cat(c("dxy cutoff:", dxy_cutoff),
    file = paste0(outdir, "/analyses/dxy/",pop1, "_", pop2, "/", pop1, "_", pop2, "dxy_stats.txt"),
    sep = "\n",
    append = TRUE)

# Save outlier file
outlier_file <- file.path(outdir, "analyses/dxy", paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".chrom.dxy.snps.outlierdxy.csv"))
write.csv(outlier_dxy, outlier_file, row.names = FALSE)

# Prepare data for plotting
dxy$chromo <- factor(dxy$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- dxy %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(dxy, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
ggplot(plot_data, aes(x = BPcum, y = dxy)) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length(unique(dxy$chromo)) / 2)) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "Chromosome", y = "dxy") +
  geom_hline(yintercept = dxy_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = file.path(outdir, "analyses/dxy", paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".dxy.snps.sigline.png")), 
       width = 20, height = 5, units = "in")
