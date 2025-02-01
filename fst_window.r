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

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
win <- args[2]
pop1 <- args[3]
pop2 <- args[4]
color1 <- args[5]
color2 <- args[6]
cutoff <- args[7]

fst_file <- file.path(outdir, "analyses/fst", paste0(win, "/", pop1, "_", pop2, "/slidingwindow.", pop1, "_", pop2, ".chroms.txt"))
fst <- read.csv(fst_file, sep = '\t') %>% na.omit()

min_fst <- min(fst$fst)
max_fst <- max(fst$fst)

cat(c("Min fst:", min_fst),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt"),
    sep = "\n", append = TRUE)

cat(c("Max fst:", max_fst),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt"),
    sep = "\n", append = TRUE)


# z transform fst values

fst_xbar <- mean(fst$fst, na.rm = TRUE)
fst_sd <- sd(fst$fst, na.rm = TRUE)
fst$z <- (fst$fst - fst_xbar) / fst_sd
p_values_one_tailed <- pnorm(q = fst$z, lower.tail = FALSE)

# Calculate the -log10 of the p-value
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

# Identify top outliers
nsnps <- nrow(fst)
top_snps <- round(nsnps * cutoff)
outlier_fst <- fst %>% 
  arrange(desc(neg_log_pvalues_one_tailed)) %>% 
  head(top_snps) %>% 
  arrange(chr, midPos)


# Save cutoff value
cat(c("fst cutoff:", fst_cutoff),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt"),
    sep = "\n",
    append = TRUE)

# Save outlier file
outlier_file <- file.path(outdir, "analyses/fst", paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".chrom.fst.snps.outlierfst.csv"))
write.csv(outlier_fst, outlier_file, row.names = FALSE)

# Prepare data for plotting
fst$chr <- factor(fst$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- fst %>%
  group_by(chr) %>%
  summarise(chr_len = max(midPos)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(fst, by = "chr") %>%
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

axisdf <- plot_data %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum))

# Plot
ggplot(plot_data, aes(x = BPcum, y = fst)) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length(unique(fst$chr)) / 2)) +
  scale_x_continuous(labels = axisdf$chr, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "Chromosome", y = "fst") +
  geom_hline(yintercept = fst_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = file.path(outdir, "analyses/fst", paste0(pop1, "_", pop2, "/", pop1, "_", pop2, ".fst.snps.sigline.png")), 
       width = 20, height = 5, units = "in")
