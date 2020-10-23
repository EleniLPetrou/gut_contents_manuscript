# The purpose of this script is to plot the results of rubias GSI simulations
# Data analyzed using "R version 3.6.1 (2019-07-05)"

################################################################################
# Load libraries

# Load libraries
library(tidyverse)
library(cowplot)

################################################################################
# Specify the input file names
hundysim_3repgrp_file <- "./analysis_rubias/analyses_manuscript_simulations/results_rubias_analysis_self-ass_100sim_3repgroups.txt"
mixsim_3repgrp_file <- "./analysis_rubias/analyses_manuscript_simulations/results_rubias_analysis_self-ass_mixsim_3repgroups.txt"
validation_file <- "./analysis_rubias/analyses_manuscript_baseline_validation/results_rubias_analysis_baseline_validation_3repgroups.txt"

################################################################################
# Read in data
hundysim_3repgrp_df <- read.delim(hundysim_3repgrp_file, header = TRUE, sep = "\t")
mixsim_3repgrp_df <- read.delim(mixsim_3repgrp_file, header = TRUE, sep = "\t")
val_df<- read.delim(validation_file, header = TRUE, sep = "\t")

################################################################################
# Plot the data
my_cols <- c("#2b83ba","#abdda4",  "#fdae61") #set some colors
mylabel1 = c(Jan_Feb = "Jan-Feb", Mar_Apr=  "Mar-Apr", May_Jun = "May")

# Make a plot of the GSI simulation results - linear regression
plotA <- ggplot(mixsim_3repgrp_df, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Simulated mixing proportion") +
  ylab("Estimated mixing proportion")+
  facet_wrap(~ repunit, labeller=labeller(repunit = mylabel1)) +
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_colour_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols)

plotA

mixsim_3repgrp_df$repunit


# Make a boxplot of GSI 100% simulations

plotB <- ggplot(hundysim_3repgrp_df, aes(x = repunit_scenario, y = reprop_posterior_mean))+
  geom_boxplot(aes(color = repunit))+
  ylim(0.5,1) +
  theme_bw()+
  xlab("Simulated proportion: 100%") +
  ylab("Estimated mixing proportion")+
  scale_x_discrete(labels = c("Jan-Feb", "Mar-Apr", "May"))+
  theme(legend.position = "none")+
  scale_color_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols)

plotB


# Make unstacked barplots with credible intervals, using the results of empirical validation of baseline

# Set some labels for plot
mylabels = c(ElBy13 = "Elliot Bay, April 2013", Marr15=  "Cherry Point, May 2014", Squa14 = "Squaxin Pass, February 2014")

# Order the label factor in a specific way

val_df$mixture_collection <- factor(val_df$mixture_collection, levels = c("Squa14", "ElBy13", "Marr15"))

# unstacked barplots, with error bars

plotC <- ggplot(data = val_df, 
                          aes(x = "", 
                              y = repprop,
                              ymin= loCI, 
                              ymax= hiCI,
                              fill = repunit)) +
  geom_bar( stat = "identity", alpha = 0.8, position = "dodge")+
  geom_errorbar(position = position_dodge(), colour="dark grey") +
  ylab("Estimated mixing proportion") +
  xlab ("") +
  facet_wrap(~mixture_collection, labeller=labeller(mixture_collection = mylabels))+
  theme_classic()+
  scale_fill_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols)

plotC


# Make a multiplot and include assessment plot with empirical samples!!
manuscript_plot <- ggdraw() +
  draw_plot(plotA, x = 0, y = 0.5, width = 0.50, height = 0.45) +
  draw_plot(plotB, x = 0.55, y = 0.5, width = 0.40, height = 0.45)+
  draw_plot(plotC, x = 0.0, y = 0.0, width = 1, height = 0.45)+
  draw_plot_label(label = c("A","B", "C"), size = 14, x = c(0, 0.55, 0), y = c(1,1,0.5))

manuscript_plot

# Save the plot to a pdf file
ggsave("./analysis_rubias/multiplot_GSI_simulations_and_validation.pdf", manuscript_plot)
