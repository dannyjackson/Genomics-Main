#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

outdir <- args[1]
win <- args[2]
pop1 <- args[3]
pop2 <- args[4]

library(qqman)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


fst <-
  read.csv(paste0(outdir, "/analyses/fst/singlesnps_fst.",
                  pop1, "_", pop2, ".txt"),
           sep = "\t")

fst_no_na <- na.omit(fst)
nrow(fst) - nrow(fst_no_na)
fst <- fst_no_na
min_fst <- min(fst$fst)
max_fst <- max(fst$fst) # print to file
cat(c("Min FST cutoff:", min_fst),
    file = "FST_stats.txt",
    sep = "\n", append = TRUE)

cat(c("Max FST cutoff:", max_fst),
    file = "FST_stats.txt",
    sep = "\n", append = TRUE)

# z transform fst values

fst_xbar <- mean(fst$fst)
fst_sd <- sd(fst$fst)

fst$z <- (fst$fst - fst_xbar) / fst_sd

p_values_one_tailed <- pnorm(q = fst$z, lower.tail = FALSE)

# Calculate -log10 of the p-value
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

ordered_fst <- fst %>%
  # desc orders from largest to smallest
  arrange(desc(neg_log_pvalues_one_tailed))

nsnps <- nrow(ordered_fst)
top_snps <- round(nsnps * 0.05)

outlier_fst_disorder <- ordered_fst[1:top_snps, ]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

fst_cutoff <- min(outlier_fst_disorder$fst) # print to file

cat(c("FST cutoff:", fst_cutoff),
    file = "FST_stats.txt",
    sep = "\n",
    append = TRUE)

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select <- -c(region))

write.csv(outlier_fst_disorder2,
          paste0(outdir, "analyses/fst/",
                 pop1, "_", pop2, ".outlierfst.csv"))

# draw it with cutoff line

blues <- c("#4EAFAF", "#082B64")

middlechr <- (max(fst$midPos) + as.numeric(win) / 2) / 2

png(file <- paste0(outdir, "/analyses/fst/",
                   pop1, "_", pop2, ".fst.snps.sigline.png"),
    width <- 2000, height = 500)

ggplot(fst, aes(x = midPos, y = (fst))) +
  # Show all points
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 0.5) +
  scale_color_manual(values <- blues) +

  # custom X axis:
  # expand=c(0, 0) removes space between plot area and x axis
  scale_y_continuous(expand <- c(0, 0), limits <- c(0, 1)) +

  # add plot and axis titles
  ggtitle(NULL) +
  labs(x <- "Chromosome", y <- "fst") +

  # add genome-wide sig and sugg lines
  geom_hline(yintercept <- fst_cutoff) +

  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data = fst[fst$is_annotate == "yes", ],
                   aes(label = as.factor(midPos), alpha = 0.7),
                   size = 5, force = 1.3) +

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