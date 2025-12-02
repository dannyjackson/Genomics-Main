library(dplyr)
library(ggplot2)

out_path_pre = '~/Desktop/darwin_finches/gone/Output_Ne_cra_pre_1'
out_path_post = '~/Desktop/darwin_finches/gone/Output_Ne_cra_post_2'

df_pre <- read.csv(out_path_pre, sep = '\t')
df_post <- read.csv(out_path_post, sep = '\t')

df_pre <- df_pre %>% mutate(time='pre')

df_post <- df_post %>% mutate(time='post')

df_full <- rbind(df_pre, df_post)
df_full <- df_full %>% 
  filter(Generation <=200)

df_full %>% 
  ggplot(aes(x=Generation, y=Geometric_mean, fill = time)) + 
  geom_point() #+coord_cartesian(ylim = c(0,20)) 
