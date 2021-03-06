# This R  script was written by Eleni Petrou on 20191003 using  "R version 3.6.1 (2019-07-05)"
# Its purpose is to make a map of sampling locations using GPS data

###########################################################
# Load R libraries
library(tidyverse)
library(maps)
library(maptools)
library(viridis)
library(ggplot2)
library(grid)
library(cowplot)

###########################################################
# Specify the file name for the input file
fileName <- "./data_sample.metadata/metadata_samples_LLTK1toLLTK8_SalishSea.txt"

##########################################################
#Read in the data
gut_df <- read.delim(fileName, header = TRUE, sep = "\t")
head(gut_df)

##########################################################
#Process the data for plotting

# Remove NAs from the data for plotting and
# Subset the data frame so you are plotting a unique Salmon ID each time
filt_df <- gut_df %>% 
  drop_na(month, latitude, longitude, julian_date, maturity)%>%
  group_by(salmon_ID) %>%
  slice(1)

head(filt_df)


# Subset the data such that you separate out the two size classes of salmon

juvs <- filt_df %>%
  filter(maturity == "juvenile")

adults <- filt_df %>%
  filter(maturity == "adult")

#########################################################
# plot the data
# Get the world polygon and extract USA and Canada

USA <- map_data("world") %>% 
  filter(region=="USA")

Canada <- map_data("world") %>% 
  filter(region=="Canada")

Mexico<- map_data("world")%>% 
  filter(region=="Mexico")

# Set the breaks for your color ramp and define labels
mybreaks=c(0,30, 60,90, 120,150, 180,210, 240,270, 300)

mylabels=c("January", "February", "March",
           "April", "May","June", 
           "July", "August", 
           "September",  "October", "November")

# Plot all samples together
plot_all <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data= filt_df, 
             aes(x=round(longitude,1), y= round(latitude,1), color= maturity), 
             alpha = 0.7) +
  coord_map(xlim= c(-121, -126),  ylim = c(47,51))+
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  labs(x = "Longitude", y = "Latitude")+
  facet_wrap(~month) +
  scale_size_continuous(name = "Number of salmon")

plot_all


##################################################################################

# Plot juvenile salmon only

plot_juvs <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data= juvs, aes(x=round(longitude,2), y= round(latitude,2), color= julian_date), alpha = 0.8, size = 2) +
  geom_point(data= adults, aes(x=round(longitude,2), y= round(latitude,2), color= julian_date), alpha = 0, size = 2) +
  coord_map(xlim= c(-122, -125),  ylim = c(47, 49))+
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) + 
  labs(x = "Longitude", y = "Latitude") + 
  scale_color_viridis(option="viridis",
                      name="Sampling date",
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, 
                      end = 1)+
  scale_size_continuous(name = "Number of salmon")

plot_juvs


# Plot only adult salmon

plot_adults <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data= juvs, aes(x=round(longitude,2), y= round(latitude,2), color= julian_date), alpha = 0, size = 2) +
  geom_point(data= adults, aes(x=round(longitude,2), y= round(latitude,2), color= julian_date), alpha = 0.6, size = 2) +
  coord_map(xlim= c(-121, -126),  ylim = c(47,51))+
  annotate("text", x = -123, y = 50, label = "Canada", size = 4) +
  annotate("text", x = -121.5, y = 48, label = "United States", size = 4)+
  labs(x = "Longitude", y = "Latitude") + 
  scale_color_viridis(option="viridis",
                      name="Sampling date",
                      breaks = mybreaks,
                      labels = mylabels,
                      begin = 0, 
                      end = 1)+
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        legend.title = element_text(size= 13),
        legend.text = element_text(size= 11)) 
  
plot_adults


# Make a little insert plot of North America

insert_map <- ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3) +
  geom_polygon(data = Mexico, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3) +
  coord_map(xlim= c(-140, -55),  ylim = c(58,20))+
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=11),
        legend.text = element_text(size=10))+
  annotate("rect", xmin = -121, xmax = -125, ymin = 47, ymax = 50, color = "black", alpha = 0)

insert_map


# Plot with multiple panels

multiplot <- ggdraw() +
  draw_plot(plot_adults, x = 0, y = 0.1, width = 0.8, height = 0.8) + #adult plot
  draw_plot(plot_juvs + theme(legend.position="null"), x = 0.55, y = 0.62, width = .28, height = .28) +
  draw_plot(insert_map, x = 0.59, y = 0.10, width = .24, height = .25) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0.10, 0.58, 0.58), y = c(0.95, 0.95, 0.38))

multiplot


# Save plot to pdf

ggsave("map_manuscript.pdf", plot = multiplot)


