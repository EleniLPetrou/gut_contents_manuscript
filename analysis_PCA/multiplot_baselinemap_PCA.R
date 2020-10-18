# The purpose of this script is to make a multiplot combining a map of the sampling
# locations of spawning herring used in the genetic baseline with the results of a PCA

# Libraries
library(tidyverse)
library(maps)
library(maptools)
library(ggrepel)
library(viridis)
library(pals)
library(grid)
library(cowplot)

##########################################################
# Specify the file names for the input files

baseline_file <- "./data_sample.metadata/metadata_baseline.txt"

PCA_baseline <-  "./analysis_PCA/results_PCA_baseline_samples.txt"

PCA_guts<- "./analysis_PCA/results_PCA_gut_samples.txt"


##########################################################
# Read in the data

baseline_df <- read.delim(baseline_file, header = TRUE, sep = "\t", comment.char = "#")

baseline_PCA_df <- read.delim(PCA_baseline, header = TRUE, sep = "\t", comment.char = "#")

gut_PCA_df <- read.delim(PCA_guts, header = TRUE, sep = "\t", comment.char = "#")


##########################################################
# Process the data
# Remove temporally replicated samples collected at the same location, so there is no
# overplotting on the map

baseline_df <- baseline_df %>%
  filter(population != "Squa14" & population != "ChPt14")

##########################################################
# Get the world polygon and extract USA and Canada

USA <- map_data("world") %>% 
  filter(region=="USA")

Canada <- map_data("world") %>% 
  filter(region=="Canada")

Mexico <- map_data("world") %>% 
  filter(region=="Mexico")


##########################################################
# Plot the data

# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150)
mylabels = c("January", "February", "March", "April", "May", "June")

# Genetic baseline map
plotA <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_label_repel( data= baseline_df, aes(x=longitude, y=latitude, label= code), size= 4) +
  geom_point(data= baseline_df, aes(x=longitude, y=latitude, color= julian_date), size = 5, alpha = 1) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size= 12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  coord_map(xlim= c(-121, -125),  ylim = c(47,50)) +
  labs(x = "Longitude", y = "Latitude") +
  
  scale_color_viridis(option="plasma",
                      name="Sampling date",
                      breaks = mybreaks, 
                      labels = mylabels,
                      begin = 0, end = 0.9) +
  panel_border(colour = "black", size = 0.5, linetype = 1,
               remove = FALSE)+
  theme(legend.position = "none")


plotA

# PCA

plotB <- ggplot() +
  geom_point(data = baseline_PCA_df, aes(x = Axis1, y = Axis2, color = julian_date), 
             size = 4, alpha = 0.4, shape = 16)+
  geom_point(data = gut_PCA_df, aes(x = Axis1, y = Axis2, shape = population), 
             size = 3, alpha = 0.6)+
  scale_color_viridis(option="plasma", name="Collection date of \nspawning herring", labels = mylabels,
                      breaks = mybreaks, begin = 0, end = 0.9)+
  scale_shape_manual(values = c(1, 2), name = "Geographic origin of\ngut content sample", labels = c("USA", "Canada"))+
  #theme(axis.title = element_blank())+
  theme_bw()+
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 11))+
  ylim(-6,3)+
  xlab("30% of variation")+
  ylab("15% of variation")

plotB

legend <- get_legend(plotB) #save the legend as an object

plotB <- plotB+ theme(legend.position = "none")

# combine the plots
final_plot<- ggdraw() +
  draw_plot(plotA, x = -0.15, y = 0.5, width = .75, height = .5)+
  draw_plot(legend, x = 0.20, y = 0.52, width = .75, height = .5)+
  draw_plot(plotB, x = 0.03, y = 0, width = .6, height = .5) +
  draw_plot_label(label = c("A", "B"), size = 14, x = c(0, 0), y = c(1, 0.52))

final_plot

##########################################################
# Save the plot to a pdf file
ggsave("./analysis_PCA/multiplot_baseline_map_and_PCA.pdf", plot = final_plot)
