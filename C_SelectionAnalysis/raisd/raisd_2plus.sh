





# RAiSD Output

# noca urban
cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/urban

ls RAiSD_Report* | awk 'BEGIN {FS = "."} {print $3"."$4}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt

cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/RAiSD_Report.noca_urban.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/urban/RAiSD_Report.noca_urban."$chrom" >> /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/RAiSD_Report.noca_urban.chromosomes

done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt





cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/

cp RAiSD_Report.noca_urban.chromosomes RAiSD_Report.noca_urban.chromosomes.plot

sed -i 's/VYXE//g' RAiSD_Report.noca_urban.chromosomes.plot
sed -i 's/\.1 /\t/g' RAiSD_Report.noca_urban.chromosomes.plot

R
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.noca_urban.chromosomes.plot', sep ='\t')
raisd_noNA <- na.omit(raisd)
nrow(raisd) - nrow(raisd_noNA)

max(raisd$U)
# 26.33
min(raisd$U)
# 3.041e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 1100740 snps, so top 0.1% would be 1100

outlier_u_disorder <- ordered_u[1:1100,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 6.253
max(outlier_u_disorder$U)
# 26.33
write.csv(outlier_u, "noca_urban.outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.noca_urban.chromosomes.plot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("noca_urban.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 2038 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 6.253) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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


# noca rural

cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/rural

ls RAiSD_Report* | awk 'BEGIN {FS = "."} {print $3"."$4}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt

cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/RAiSD_Report.noca_rural.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/rural/RAiSD_Report.noca_rural."$chrom" >> /xdisk/mcnew/dannyjackson/cardinals/raisd/noca/RAiSD_Report.noca_rural.chromosomes

done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt

cp RAiSD_Report.noca_rural.chromosomes RAiSD_Report.noca_rural.chromosomes.plot

sed -i 's/VYXE//g' RAiSD_Report.noca_rural.chromosomes.plot
sed -i 's/\.1 /\t/g' RAiSD_Report.noca_rural.chromosomes.plot

module load R/4.4.0

R
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.noca_rural.chromosomes.plot', sep ='\t')
raisd_noNA <- na.omit(raisd)
nrow(raisd) - nrow(raisd_noNA)

max(raisd$U)
# 23.84
min(raisd$U)
# 1.138e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 1023066 snps, so top 0.1% would be 1023

outlier_u_disorder <- ordered_u[1:1023,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 6.509
max(outlier_u_disorder$U)
# 23.84
write.csv(outlier_u, "noca_rural.outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.noca_rural.chromosomes.plot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("noca_rural.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 2038 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 6.509) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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
















# pyrr urban
cd /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/urban

ls RAiSD_Report* | awk 'BEGIN {FS = "."} {print $3"."$4}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt

cd /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/RAiSD_Report.pyrr_urban.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/urban/RAiSD_Report.pyrr_urban."$chrom" >> /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/RAiSD_Report.pyrr_urban.chromosomes

done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt



cp RAiSD_Report.pyrr_urban.chromosomes RAiSD_Report.pyrr_urban.chromosomes.plot

sed -i 's/VYXE//g' RAiSD_Report.pyrr_urban.chromosomes.plot
sed -i 's/\.1 /\t/g' RAiSD_Report.pyrr_urban.chromosomes.plot

R
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.pyrr_urban.chromosomes.plot', sep ='\t')
raisd_noNA <- na.omit(raisd)
nrow(raisd) - nrow(raisd_noNA)

max(raisd$U)
# 51.49
min(raisd$U)
# 4.494e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 1031424 snps, so top 0.1% would be 1031

outlier_u_disorder <- ordered_u[1:1031,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 6.796
max(outlier_u_disorder$U)
# 51.49
write.csv(outlier_u, "pyrr_urban.outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.pyrr_urban.chromosomes.plot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("pyrr_urban.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 2038 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 6.796) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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










# pyrr rural
cd /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/rural

ls RAiSD_Report* | awk 'BEGIN {FS = "."} {print $3"."$4}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt

cd /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/RAiSD_Report.pyrr_rural.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/rural/RAiSD_Report.pyrr_rural."$chrom" >> /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr/RAiSD_Report.pyrr_rural.chromosomes

done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/chromlist.txt


cp RAiSD_Report.pyrr_rural.chromosomes RAiSD_Report.pyrr_rural.chromosomes.plot

sed -i 's/VYXE//g' RAiSD_Report.pyrr_rural.chromosomes.plot
sed -i 's/\.1 /\t/g' RAiSD_Report.pyrr_rural.chromosomes.plot

R
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.pyrr_rural.chromosomes.plot', sep ='\t')
raisd_noNA <- na.omit(raisd)
nrow(raisd) - nrow(raisd_noNA)

max(raisd$U)
# 31.18
min(raisd$U)
# 4.028e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 1012748 snps, so top 0.1% would be 1013

outlier_u_disorder <- ordered_u[1:1100,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 6.522
max(outlier_u_disorder$U)
# 31.18
write.csv(outlier_u, "pyrr_rural.outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.pyrr_rural.chromosomes.plot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("pyrr_rural.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 2038 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 6.522) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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



# generate relevant gene lists
# noca
cd /xdisk/mcnew/dannyjackson/cardinals/raisd/noca

awk 'BEGIN {FS = ","} {$1=""}1' noca_rural.outlieru.tsv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > noca_rural.outlieru.tmp.tsv

mv noca_rural.outlieru.tmp.tsv noca_rural.outlieru.tsv

awk 'BEGIN {FS = ","} {$1=""}1' noca_urban.outlieru.tsv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > noca_urban.outlieru.tmp.tsv

mv noca_urban.outlieru.tmp.tsv noca_urban.outlieru.tsv

sed -i 's/\"//g' noca_urban.outlieru.tsv
sed -i 's/\"//g' noca_rural.outlieru.tsv

# pyrr
cd /xdisk/mcnew/dannyjackson/cardinals/raisd/pyrr
awk 'BEGIN {FS = ","} {$1=""}1' pyrr_rural.outlieru.tsv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > pyrr_rural.outlieru.tmp.tsv

mv pyrr_rural.outlieru.tmp.tsv pyrr_rural.outlieru.tsv

awk 'BEGIN {FS = ","} {$1=""}1' pyrr_urban.outlieru.tsv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > pyrr_urban.outlieru.tmp.tsv

mv pyrr_urban.outlieru.tmp.tsv pyrr_urban.outlieru.tsv

sed -i 's/\"//g' pyrr_rural.outlieru.tsv
sed -i 's/\"//g' pyrr_urban.outlieru.tsv


# noca urban
while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_noca_urban_raisd.txt

done < noca_urban.outlieru.tsv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_urban_raisd.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_urban_raisd.txt


# noca rural

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_noca_rural_raisd.txt

done < noca_rural.outlieru.tsv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_rural_raisd.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_rural_raisd.txt



# pyrr urban
while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_pyrr_urban_raisd.txt

done < pyrr_urban.outlieru.tsv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_urban_raisd.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_urban_raisd.txt


# pyrr rural

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_pyrr_rural_raisd.txt

done < pyrr_rural.outlieru.tsv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_rural_raisd.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_rural_raisd.txt
