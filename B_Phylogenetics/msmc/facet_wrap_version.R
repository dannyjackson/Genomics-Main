library(ggplot2)
library(dplyr)
library(tidyverse)
library(forcats)
library(patchwork)
library(ggh4x)



mu = 1.0e-8
gen = 4.5
out_path = '~/Desktop/new_msmc_outputs/lines/all'

#for (spp_folder in list.files(out_path)) {

# Get output file list for species
spp_path <- paste0(out_path, spp_folder, '/')
files <- rev(list.files(path = out_path, full.names = TRUE))

# Set future variable name for plot object
plot_var_name <- paste0(spp_folder, '_plot')

# Setting some plot parameters we want depending on spp
#if (spp_folder == 'cra') {
#  post_color = "#4EAFAF"
#  min_x_lim = 1200
#  plot_label = 'A'
#} else if (spp_folder == 'par') {
#  post_color = "#A6C965"
#  min_x_lim = 1200
#  plot_label = 'B'
#} else if (spp_folder == 'for') {
#  post_color = "#FF817E"
#  min_x_lim = 4000
#  plot_label = 'C'
#}


df_lst <- lapply(files, function(f) {
  df <- read.csv(f, sep = '\t')
  df$id <- basename(f)   # add file name as id
  df$period <- strsplit(basename(f), "[_.]+")[[1]][4] # color based on whether file is pre or post. Note that this string parsing may change depending on your filenames
  df$type <- strsplit(basename(f), "[_.]+")[[1]][5]# Make column stating if entry is from data or bootstraps. Makes plotting varying alpha values easier
  df$spp <- strsplit(basename(f), "[_.]+")[[1]][3]# Make column stating spp. Separating out data for facet_wrap
  df$spp_period <- paste0(strsplit(basename(f), "[_.]+")[[1]][3], '_', strsplit(basename(f), "[_.]+")[[1]][4])
  return(df)
})

combined_df <- do.call(rbind, df_lst)

combined_df <- combined_df %>% 
  mutate(ne = (1/lambda)/(2*mu),
         time_years = (left_time_boundary * gen) / mu,
         ne_thousands = ((1/lambda)/(2*mu)) /1000 )

# Do this to ensure that the post population (the one we're most interested in is plotted in front)
combined_df$id <- fct_inorder(combined_df$id)

color_vect <- c('PAR_pre'='gray41', 'FOR_pre'='gray41','CRA_pre'='gray41',
                'PAR_post'="#A6C965", 'FOR_post'="#FF817E", 'CRA_post'="#4EAFAF")

custom_x_scales <- list(
  scale_x_continuous(limits = c(1100, 1000000)),
  scale_x_continuous(limits = c(4000, 1000000)),
  scale_x_continuous(limits = c(1100, 1000000))
)


p = ggplot(combined_df,aes(y = ne_thousands,x = time_years, group=id, color=spp_period, alpha = type, linewidth = type)) + 
  geom_step(key_glyph='rect') +
  facet_wrap(~ spp, nrow=3, scales='free') +
  scale_color_manual(labels=c('Post-Invasion','Pre-Invasion'), values = color_vect) +
  scale_alpha_manual(values = c('data'=1, 'bootstrap'=0.4)) +
  scale_linewidth_manual(values = c('data'=2, 'bootstrap'=0.6)) +
  coord_cartesian(ylim = c(0,150), xlim = c(min_x_lim, 1000000)) + 
  scale_x_log10() + 
  labs(x="Years Before Present (log scale)", y="Effective Population Size (in thousands of individuals)", color='', title = plot_label) +
  guides(alpha='none', linewidth='none',
         color = guide_legend(override.aes = list(15))) +
  theme(axis.title = element_text(size = 30), plot.title = element_text(size = 35),
        legend.text = element_text(size = 25), axis.text = element_text(size = 25), 
        legend.position = c(0.18,0.9), legend.key.size = unit(1.6,'cm'),
        plot.margin = unit(c(0,20,10,20), 'pt'), axis.title.x = element_text(vjust = -0.5),
        panel.spacing.y = unit(2, "lines")) +
  facetted_pos_scales(x=custom_x_scales)

#assign(plot_var_name, p)

#}

#all_msmc_plots <- cra_plot + for_plot + par_plot
print(p)
ggsave('~/Desktop/new_msmc_outputs/msmc_plots.png', plot=p, , width=10, height=18, units='in')



