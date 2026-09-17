# Basic example of plotting gone results across multiple populations in a time series
library(dplyr)
library(ggplot2)

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
out_path_pre <- args[1] # Usually the earlier ("pre") population
out_path_post <- args[2] # Usually the later ("post") population
num_gens <- args[3] # Number of generations
color_pre <- args[4]
color_post <- args[5]

 
df_pre <- read.csv(out_path_pre, sep = '\t')
df_post <- read.csv(out_path_post, sep = '\t')

df_pre <- df_pre %>% mutate(time='pre') # Add time info
df_post <- df_post %>% mutate(time='post')

df_full <- rbind(df_pre, df_post) # Combine data comparing Ne estimates between time-separated populations
df_full <- df_full %>%
  filter(Generation <=num_gens) # Apply whatever data filters you want. GONE is only good up to 200 generations

color_codes <- c(pre=color_pre, post=color_post)

df_full %>%
  ggplot(aes(x=Generation, y=Ne_diploids, color = time)) +
  geom_step() +
  scale_color_manual(values = color_codes, name='Time') #+coord_cartesian(ylim = c(0,20))
