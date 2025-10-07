library(ggplot2)
library(dplyr)
library(tidyverse)
library(forcats)
library(patchwork)

# This script is designed where your msmc output files for each species are present in their own directory underneath a single "out_path" directory
# In a loop, we will plot msmc outputs for each species directory

# SCRIPT PARAMETERS ===========================================
mu = 1.0e-8
gen = 4.5
out_path = 'data/'

#pdf('2025-09-04_msmc_plots.pdf', width = 15, height = 20)

# =============================================================

for (spp_folder in list.files(out_path)) {

  # Get output file list for species
  spp_path <- paste0(out_path, spp_folder, '/')
  print(spp_folder)
  files <- rev(list.files(path = spp_path, full.names = TRUE))

  # Set future variable name for plot object
  plot_var_name <- paste0(spp_folder, '_plot')

  # Setting some plot parameters we want depending on spp
  if (spp_folder == 'cra') {
    post_color = "#4EAFAF"
    min_x_lim = 7000
    plot_label = 'A'
  } else if (spp_folder == 'for') {
    post_color = "#FF817E"
    min_x_lim = 7000
    plot_label = 'B'
  } else if (spp_folder == 'par') {
    post_color = "#A6C965"
    min_x_lim = 7000
    plot_label = 'C'
  }


  df_lst <- lapply(files, function(f) {
    df <- read.csv(f, sep = '\t')
    df$id <- basename(f)   # add file name as id
    df$period <- strsplit(basename(f), "[_.]+")[[1]][4] # color based on whether file is pre or post. Note that this string parsing may change depending on your filenames
    df$type <- strsplit(basename(f), "[_.]+")[[1]][5]# Make column stating if entry is from data or bootstraps. Makes plotting varying alpha values easier
    return(df)
  })

  #Stack all dataframes and insert needed data columns
  combined_df <- do.call(rbind, df_lst)
  combined_df <- combined_df %>%
    mutate(ne = (1/lambda)/(2*mu),
           time_years = (left_time_boundary * gen) / mu,
           ne_thousands = ((1/lambda)/(2*mu)) /1000 ) # Note that this column is the same as ne, but stores ne values in thousands

  # Do this to ensure that the post population (the one we're most interested in) is plotted in front
  combined_df$id <- fct_inorder(combined_df$id)

  color_vect <- c('pre'='gray41', 'post'=post_color)


  p = ggplot(combined_df,aes(y = ne_thousands,x = time_years, group=id, color=period, alpha = type, linewidth = type)) +
    geom_step(key_glyph='rect') + # key_glyph sets legend symbol shape
    # Set custom line parameters
    scale_color_manual(labels=c('Post-Invasion','Pre-Invasion'), values = color_vect) +
    scale_alpha_manual(values = c('data'=1, 'bootstrap'=0.4)) +
    scale_linewidth_manual(values = c('data'=2, 'bootstrap'=0.6)) +
    # Set scale limits/parameters
    coord_cartesian(ylim = c(0,150), xlim = c(min_x_lim, 1000000)) +
    scale_x_log10() +
    labs(x="Years Before Present (log scale)", y=bquote("Effective Population Size" ~(x10^3)), color='', title = plot_label) +
    # Do this to remove extra unneeded legends
    guides(alpha='none', linewidth='none',
           color = guide_legend(override.aes = list(15))) +
    # Lots of final visual touch-ups
    theme_bw() +
    theme(axis.title = element_text(size = 22), plot.title = element_text(size = 35), legend.text = element_text(size = 15),
          axis.text = element_text(size = 20), legend.position = c(0.18,0.85), legend.key.size = unit(1,'cm'),
          plot.margin = unit(c(10,20,10,20), 'pt'), axis.title.x = element_text(vjust = -0.5))

  assign(plot_var_name, p)

}

all_msmc_plots <- cra_plot + for_plot + par_plot

ggsave('2025-09-04_msmc_plots_square.pdf', plot=all_msmc_plots, , width=25, height=8, units='in', device = 'pdf')
#dev.off()
