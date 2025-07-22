#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)
#OUTDIR = args[1]
#IND = args[2]

# Add some code here based on your needs to both parse your directories for the msmc outputs you need and generate the types of plots you want
#For example, the lazy way I've done it is shown below

# This script only works currently on a local machine in RStudio. Not yet set up for submitting via sbatch on HPC. May implement this later if feel like it.
# Generating the plots is simple and not resource-intensive, and you'll likely want to adjust your plots frequently (which is easier in live environment like RStudio)


# All post-philornis samples
post_lst = c('JP4481_all','JP5410_all','JP9655_all','SM1067_all', 'SM1157_all',
             'SM1200_all','SM1231','SM1240_all','SM1266_all','SM1083','SM1156','SM1204_all','SM1237',
             'SM1270','SM1271','SM1272','SM1273','RHC097_all','RHC507_all',
             'SM031_all','SM032_all','SM040_all','SM059_all','SM079_all')
# All CRA Individuals
cra_lst = c('JP4481_all','JP5410_all','JP9655_all','lamich_PL15','lamich_PL16','lamich_PL4','lamich_PL7','lamich_PL9','SM1067_all','SM1157_all',
            'SM1200_all','SM1231','SM1240_all','SM1266_all')
# All FOR Individuals
for_lst = c('SM1083','SM1156','SM1204_all, SM1237','SM1270','SM1271','SM1272','SM1273','SRR2917289','SRR2917290',
            'SRR2917291','SRR2917292','SRR2917293','SRR2917294','SRR2917295','SRR2917296','SRR2917297','SRR2917298')
# All PAR Individuals
par_lst = c('lamich_PARV1','lamich_PARV2','RHC097_all','RHC507_all','SM031_all','SM032_all','SM040_all','SM059_all','SM079_all','SRR2917329',
            'SRR2917330','SRR2917331','SRR2917332','SRR2917333','SRR2917334','SRR2917335','SRR2917336','SRR2917337','SRR2917338')

# Plot Parameters ==============================================================
mu = 2.04e-09
gen = 5
out_path = 'path/to/msmc/outputs/'
plot_title = '__________ Ne Over Time'
pre_post_plot = FALSE
plot_type = 'bootstrapped' # Options can include pre_post, spp, bootstrapped
                       # pre_post = color lines by whether individual is pre or post philornis
                       # spp = color lines by spp of individual
                       # bootstrapped = format plot to show off bootstraps (requires stating bootstrapped_indv_fname variable)
original_data_fname = 'msmc_output_some_indv.final.txt' # Filename of original individual (nonbootstrapped data)
plot_fname = 'msmc_plot'

# CODE =========================================================================
plot_exists = FALSE
for (output in list.files(path = out_path)) {
  
  print(output)
  
  if (plot_type == 'pre_post') {
    label_names = c('Pre Philornis', 'Post Philornis')
    label_colors=c('purple', 'orange')
    line_width = 1.5
    # Coloring by pre/post
    if (grepl(paste(post_lst,collapse = '|'), output)) {
      color = 'orange'
    } else {
      color = 'purple'
    }
  }
  else if (plot_type == 'spp') {
    label_names = c('CRA', 'FOR', 'PAR')
    label_colors=c('darkgoldenrod1', 'deepskyblue3', 'deeppink3')
    line_width = 1.5
    # Coloring by species
    if (grepl(paste(cra_lst,collapse = '|'), output)) {
      color = 'darkgoldenrod1'
    } else if (grepl(paste(for_lst,collapse = '|'), output)) {
      color = 'deepskyblue3'
    } else if (grepl(paste(par_lst,collapse = '|'), output)) {
      color = 'deeppink3'
    }
  }
  else if (plot_type == 'bootstrapped') {
    label_names = c('Boostraps', 'Data')
    label_colors=c('darkgoldenrod1', 'deepskyblue3')
    # Coloring by whether data is bootstrapped
    if (output == original_data_fname) {
      color = 'deepskyblue3'
      line_width = 3
    } else {
      color = 'darkgoldenrod1'
      line_width = 1
    }
  }
  
  data <- read.csv(paste(out_path, output, sep = ''), sep ='\t')
  
  time <- (data$left_time_boundary/mu*gen)
  pop.size <- (1/data$lambda)/(2*mu)
  
  if (plot_exists == FALSE) {
    plot(time, pop.size, type='s', 
         xlab="Years before present (log scale)", 
         ylab="Effective Population Size", log='x', ylim = c(0,500000), col=color, lwd=line_width)
    plot_exists = TRUE
  } else {
    lines(time, pop.size, type = 's', col=color, lwd=line_width)
  }
  }

legend(x='topleft',legend=label_names, fill=label_colors)
title(main=plot_title)

png(paste0(out_path, plot_fname, ".cropped.png"), width=500, height=500)

dev.off()

