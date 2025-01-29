#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pop1 <- args[2]
pop2 <- args[3]

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
           sep = "\t")

fst_no_na <- na.omit(fst)
nrow(fst) - nrow(fst_no_na)
fst <- fst_no_na

min_fst <- min(fst$fst)
max_fst <- max(fst$fst)

cat(c("Min FST cutoff:", min_fst),
    file = paste0(outdir, "/analyses/fst/FST_stats.txt"),
    sep = "\n", append = TRUE)

cat(c("Max FST cutoff:", max_fst),
    file = paste0(outdir, "/analyses/fst/FST_stats.txt"),
    sep = "\n", append = TRUE)

# z transform fst values

fst_xbar <- mean(fst$fst, na.rm = TRUE)
fst_sd <- sd(fst$fst, na.rm = TRUE)

fst$z <- (fst$fst - fst_xbar) / fst_sd

p_values_one_tailed <- pnorm(q = fst$z, lower.tail = FALSE)

# Calculate -log10 of the p-value
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

ordered_fst <- fst %>%
  # desc orders from largest to smallest
  arrange(desc(neg_log_pvalues_one_tailed))

nsnps <- nrow(ordered_fst)
top_snps <- round(nsnps * 0.001)

outlier_fst_disorder <- ordered_fst[1:top_snps, ]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

fst_cutoff <- min(outlier_fst_disorder$fst) # print to file

cat(c("FST cutoff:", fst_cutoff),
    file = paste0(outdir, "/analyses/fst/FST_stats.txt"),
    sep = "\n",
    append = TRUE)

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2,
          paste0(outdir, "analyses/fst/singlesnps.",
                 pop1, "_", pop2, ".outlierfst.csv"))

# draw it with cutoff line

blues <- c("#4EAFAF", "#082B64")


fst$chr <- factor(fst$chr, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

df.tmp <- fst %>%

  # Compute chromosome size
  group_by(chr) %>%
  summarise(chr_len = max(midPos)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(fst, ., by = c("chr" = "chr")) %>%

  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

half_length <- ceiling(length(unique(fst$chr)) / 2)


# get chromosome center positions for x-axis
axisdf <- df.tmp %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

png(file = paste0(outdir, "/analyses/fst/",
                   pop1, "_", pop2, ".fst.snps.sigline.png"),
    width = 2000, height = 500)

ggplot(df.tmp, aes(x = BPcum, y = (fst))) +
  # Show all points
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(blues, half_length)) +
  # custom X axis:
  # expand=c(0,0)removes space between plot area and x axis
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +

  scale_y_continuous(expand <- c(0, 0), limits <- c(0, 1)) +
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = fst_cutoff) +
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"),color="orange", size=2) +
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  # Custom the theme:
  theme_bw(base_size <- 22) +
  theme(
    plot.title <- element_text(hjust <- 0.5),
    legend.position = "none",
    panel.border <- element_blank(),
    panel.grid.major.x <- element_blank(),
    panel.grid.minor.x <- element_blank()
  )


dev.off()