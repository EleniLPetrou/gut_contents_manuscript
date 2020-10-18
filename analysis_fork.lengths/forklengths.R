# This script was written by Eleni on 20191004 using "R version 3.6.1 (2019-07-05)"
# Its purpose is to make some nice histograms summarizing the forklength metadata

# Load the necessary libraries
library (tidyverse)
library(cowplot)
library(viridis)
library(scales)

#######################################################################
# input file names

fileName <- "./data_sample.metadata/metadata_samples_LLTK1toLLTK8_SalishSea.txt"

#######################################################################
# read in the sample metadata for plotting
salish_df <- read.delim(fileName, sep = "\t", header = TRUE)
head(salish_df)


#######################################################################
# Remove samples with missing data and subset the data by life history stage and 

nas <- salish_df %>%
  filter(is.na(maturity))

salish_df <- salish_df %>% 
  mutate("herring_length_cm" = herring_length_mm/10) %>%
  drop_na(month)

adults <- salish_df %>%
  filter(maturity == "adult")

juvs <- salish_df %>%
  filter(maturity == "juvenile")



#######################################################################
# Plot distribution of all sample lengths

# Get some viridis colors for a discrete scale- add
# complexity to the visualizations
show_col(viridis_pal(option = "D", begin = 0, end = 1)(10))
month_colors <- (viridis_pal(option = "D", begin = 0, end = 1)(10))

unique(salish_df$month)

mylabels=c("January", "February", "March",
           "April", "May", "June", "July", 
           "August", "September")



plot1 <- ggplot(salish_df, aes(x = herring_length_cm, fill = as.factor(month))) + 
  geom_histogram()+
  scale_fill_manual(name = "Collection month", values = month_colors, labels = mylabels)+
  theme_classic()+
  xlab("Herring length (cm)")

plot1 

plot2 <- ggplot(salish_df, aes(x = as.numeric(salmon_length_cm), fill = as.factor(month))) + 
  geom_histogram()+
  scale_fill_manual(name = "Collection month",values = month_colors, labels = mylabels)+
  theme_classic()+
  xlab("Salmon length (cm)")

plot2



# combine plots into one pdf

final_plot <- ggdraw() +
  draw_plot(plot1 + theme(legend.position = "none"), x = 0, y = 0.5, width = 0.8, height = 0.45) +
  draw_plot(plot2, x = 0.0, y = 0.00, width = .96, height = .45)+
  draw_plot_label(label = c("A: Herring", "B: Salmon"), size = 12,
                  x = c(0, 0), y = c(1, 0.5))
final_plot


# Save plot to pdf file

ggsave("forklengths_manuscript.pdf", plot = final_plot)
