args = commandArgs(trailingOnly=TRUE)
OUTDIR ="~/IntroBioinformaticsWorkshop"
WIN = 1000

OUTDIR = args[1]
WIN = args[2]

library(qqman)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


fst <- read.csv(file.path(OUTDIR, "analyses/fst/singlesnps_fst_pyrr.txt"), sep ='\t')

fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA
minFST = min(fst$fst)
maxFST = max(fst$fst) # print to file
cat(c("Min FST cutoff:",minFST),file="FST_stats.txt", sep="\n", append=TRUE)
cat(c("Max FST cutoff:",maxFST),file="FST_stats.txt", sep="\n", append=TRUE)

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nsnps = nrow(ordered_fst)
top_snps = round(nsnps*0.05)

outlier_fst_disorder <- ordered_fst[1:top_snps,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

FSTcutoff = min(outlier_fst_disorder$fst) # print to file

cat(c("FST cutoff:",FSTcutoff),file="FST_stats.txt", sep="\n", append=TRUE)


outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, file.path(OUTDIR, "analyses/fst/pyrr.outlierfst.csv"))



# draw it with cutoff line 

blues <- "#082B64"

middlechr = (max(fst$midPos) + as.numeric(WIN)/2)/2

png(file = file.path(OUTDIR,"analyses/fst/pyrr.fst.snps.sigline.png"), width = 2000, height =500)

ggplot(fst, aes(x=midPos, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = blues) +

  # custom X axis:
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "fst") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = FSTcutoff) +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=fst[fst$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()