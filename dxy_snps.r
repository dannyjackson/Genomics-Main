#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pop1 <- args[2]
pop2 <- args[3]
color1 <- args[4]
color2 <- args[5]

# Package names
packages <- c("qqman", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))



# Plot dxy 

dxy <- read.csv(paste0(outdir, "/analyses/dxy", pop1, "_" pop2,'/Dxy_persite_nocaurban_nocarural.txt'), sep ='\t')


dxy_no_na <- na.omit(dxy)
nrow(dxy) - nrow(dxy_no_na) 
dxy <- dxy_no_na

# z transform dxy values

dxy_xbar <- mean(dxy$dxy, na.rm = TRUE)
dxy_sd <- sd(dxy$dxy, na.rm = TRUE)

dxy$z <- (dxy$dxy - dxy_xbar) / dxy_sd

p_values_one_tailed <- pnorm(q = dxy$z, lower.tail = FALSE)

# Calculate the -log10 of the p-value
dxy$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nsnps <- nrow(dxy)
top_snps <- round(nsnps * 0.001)

outlier_dxy_disorder <- ordered_dxy[1:topsnps, ]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromo, position)

dxy_cutoff <- min(outlier_dxy_disorder$dxy) # print to file


cat(c("dxy cutoff windowed:", dxy_cutoff),
    file = paste0(outdir, "/analyses/dxy/", ${POP1}, "_", ${POP2}, "dxy_stats.txt"),
    sep = "\n",
    append = TRUE)


write.csv(outlier_dxy,
          paste0(outdir, "/analyses/dxy/",
                 pop1, "_", pop2, "/", pop1, "_", pop2, ".chrom.dxy.windowed.outlierdxy.csv"))


# draw it with cutoff line 


colors <- c(color1, color2)

dxy$chromo <- factor(dxy$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))


df.tmp <- dxy %>%

  # Compute chromosome size
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(dxy, ., by = c("chromo" = "chromo")) %>%

  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

half_length <- ceiling(length(unique(dxy$chr)) / 2)

# get chromosome center positions for x-axis
axisdf <- df.tmp %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

png(file = paste0(outdir, "/analyses/dxy/",
                   pop1, "_", pop2, "/", pop1, "_", pop2, ".dxy.snps.sigline.png"),
    width = 2000, height = 500)

ggplot(df.tmp, aes(x = BPcum, y = (dxy))) +
  # Show all points
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(colors, half_length)) +
  # custom X axis:
  # expand=c(0,0)removes space between plot area and x axis
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +

  scale_y_continuous(expand <- c(0, 0), limits <- c(0, 1)) +
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = dxy_cutoff) +
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