# The purpose of this script is to plot the results of mixed stock analyses from rubias.  
# Data analyzed using "R version 3.6.1 (2019-07-05)"

###############################################################################################
# Load libraries
library(cowplot)
library(tidyverse)

###############################################################################################
# Specify the input file names (containing reference samples in rubias format)

adult_file <- "./analysis_rubias/analyses_manuscript_gutcontents/results_rubias_analysis_gutcontents_3repgroups_adultsalmon_byseason.txt"
juv_file <- "./analysis_rubias/analyses_manuscript_gutcontents/results_rubias_analysis_gutcontents_3repgroups_juvenilesalmon_byseason.txt"

###############################################################################################
# Read in data
adult_df <- read.delim(adult_file, header = TRUE, sep = "\t")
juv_df <- read.delim(juv_file, header = TRUE, sep = "\t")

###############################################################################################
# Plot the data

my_cols3 <- c("#2b83ba","#abdda4",  "#fdae61")

# specify some custom facet labels

juv_labels <- c(spring = "spring (N= 47)", summer = "summer (N= 40)" )

adult_labels <- c(spring = "spring (N= 260)", summer = "summer (N= 64)", winter = "winter (N= 68)" )


juv_plot <- ggplot() +
  geom_bar(data = juv_df, aes(x = "", y = repprop, fill = repunit),width = 0.4, stat = "identity", alpha = 0.8)+
  geom_errorbar(data = filter(juv_df, repunit == "May_Jun"),  aes(x="", ymin= loCI, ymax= hiCI), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  geom_errorbar(data = filter(juv_df, repunit == "Mar_Apr", mixture_collection != "spring"),  aes(x="", ymin= loCI, ymax= hiCI), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  geom_errorbar(data = filter(juv_df, repunit == "Jan_Feb", mixture_collection == "spring"),  aes(x="", ymin= 1-hiCI , ymax= 1-loCI ), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  ylab("Estimated mixing proportion") +
  xlab ("") +
  facet_wrap(~mixture_collection, labeller=labeller(mixture_collection = juv_labels))+
  theme_classic()+
  scale_fill_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols3)

juv_plot 



adult_plot <- ggplot() +
  geom_bar(data = adult_df, aes(x = "", y = repprop, fill = repunit),width = 0.4, stat = "identity", alpha = 0.8)+
  geom_errorbar(data = filter(adult_df, repunit == "May_Jun"),  aes(x="", ymin= loCI, ymax= hiCI), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  geom_errorbar(data = filter(adult_df, repunit == "Mar_Apr", mixture_collection != "spring"),  aes(x="", ymin= loCI, ymax= hiCI), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  geom_errorbar(data = filter(adult_df, repunit == "Jan_Feb", mixture_collection == "spring"),  aes(x="", ymin= 1-hiCI , ymax= 1-loCI ), colour="black", alpha=0.9, size= 0.5, width = 0.2)+
  ylab("Estimated mixing proportion") +
  xlab ("") +
  facet_wrap(~mixture_collection, labeller=labeller(mixture_collection = adult_labels))+
  theme_classic()+
  scale_fill_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols3)

adult_plot 


################################################################
# make multiplots

# Version #1
plot_grid(
  juv_plot, adult_plot,
  labels = c("A", "B"),
  nrow = 2,
  align="v"
)

# Version #2

final_plot <- ggdraw() +
  draw_plot(juv_plot, x = 0, y = 0.50, width = 0.75, height = 0.45) +
  draw_plot(adult_plot+ theme(legend.position = "none"), x = 0, y = 0, width = 0.75, height = 0.45)+
  draw_plot_label(label = c("A: juvenile salmon", "B: adult salmon"), size = 12, x = c(-0.01, -0.01), 
                  y = c(0.99,0.50))

final_plot


################################################################
# Save the output
ggsave("./analysis_rubias/analyses_manuscript_gutcontents/multiplot_results.pdf", plot = final_plot)




