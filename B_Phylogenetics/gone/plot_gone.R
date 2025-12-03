library(dplyr)
library(ggplot2)

out_path_pre = '~/Desktop/darwin_finches/gone/Output_Ne_for_pre4'
out_path_post = '~/Desktop/darwin_finches/gone/Output_Ne_for_post4'

df_pre <- read.csv(out_path_pre, sep = '\t', skip = 1) # Skip first line of output to make parsable df
df_post <- read.csv(out_path_post, sep = '\t', skip = 1)

df_pre <- df_pre %>% mutate(time='pre') # Add time info
df_post <- df_post %>% mutate(time='post')

df_full <- rbind(df_pre, df_post) # Combine for data for comparing Ne estimates between time-separated populations
df_full <- df_full %>% 
  filter(Generation <=200) # Apply whatever data filters you want. Usually GONE is only good up to 200 generations

color_codes <- c(pre='black', post='red')

df_full %>% 
  ggplot(aes(x=Generation, y=Geometric_mean, color = time)) + 
  geom_point() +
  scale_color_manual(values = color_codes, name='Time') #+coord_cartesian(ylim = c(0,20)) 
