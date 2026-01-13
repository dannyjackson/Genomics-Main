library(dplyr)
library(ggplot2)

out_path_pre = 'pre_pop_ne_output_path'
out_path_post = 'post_pop_ne_ouput_path'

df_pre <- read.csv(out_path_pre, sep = '\t')
df_post <- read.csv(out_path_post, sep = '\t')

df_pre <- df_pre %>% mutate(time='pre') # Add time info
df_post <- df_post %>% mutate(time='post')

df_full <- rbind(df_pre, df_post) # Combine for data for comparing Ne estimates between time-separated populations
df_full <- df_full %>%
  filter(Generation <=200) # Apply whatever data filters you want. GONE is only good up to 200 generations

color_codes <- c(pre='black', post='red')

df_full %>%
  ggplot(aes(x=Generation, y=Ne_diploids, color = time)) +
  geom_point() +
  scale_color_manual(values = color_codes, name='Time') #+coord_cartesian(ylim = c(0,20))

