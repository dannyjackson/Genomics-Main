#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pop1 <- args[2]
pop2 <- args[3]
color1 <- args[4]
color2 <- args[5]
cutoff <- args[6]

# Package names
packages <- c("qqman", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

fst <-
  read.csv(paste0(outdir, "/analyses/fst/singlesnps.",
                  pop1, "_", pop2, ".chroms.txt"),
           sep = "\t") %>% na.omit()

min_fst <- min(fst$fst)
max_fst <- max(fst$fst)

cat(c("Min FST:", min_fst),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt")),
    sep = "\n", append = TRUE)

cat(c("Max FST:", max_fst),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt")),
    sep = "\n", append = TRUE)

# z transform fst values
fst_xbar <- mean(fst$fst, na.rm = TRUE)
fst_sd <- sd(fst$fst, na.rm = TRUE)
fst$z <- (fst$fst - fst_xbar) / fst_sd
p_values_one_tailed <- pnorm(q = fst$z, lower.tail = FALSE)

# Calculate -log10 of the p-value
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

# Identify top outliers
nsnps <- nrow(ordered_fst)
top_snps <- round(nsnps * cutoff)
outlier_fst <- fst %>% 
  arrange(desc(neg_log_pvalues_one_tailed)) %>% 
  head(top_snps) %>% 
  arrange(chromo, position)

fst_cutoff <- min(outlier_fst_disorder$fst) # print to file

# save cutoff value
cat(c("FST cutoff:", fst_cutoff),
    file = paste0(outdir, "/analyses/fst/",pop1, "_", pop2, "/", pop1, "_", pop2, "fst_stats.txt")),
    sep = "\n",
    append = TRUE)

# Save outlier file
write.csv(outlier_fst,
          paste0(outdir, "analyses/fst/singlesnps.",
                 pop1, "_", pop2, ".outlierfst.csv"))


# Prepare data for plotting
fst$chromo <- factor(fst$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- fst %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(fst, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
ggplot(plot_data, aes(x = BPcum, y = fst)) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length(unique(fst$chromo)) / 2)) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
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
dev.off()