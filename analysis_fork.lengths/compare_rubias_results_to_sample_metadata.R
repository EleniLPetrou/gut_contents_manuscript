# This script was written by Eleni on 20200214
# It takes individual posterior probabilities from rubias and explores 
# their relationship to some basic sample metadata

##############################################################################
# Load libraries

library (tidyverse)
library(cowplot)
library(viridis)
library(scales)

##############################################################################
# Specify file names

metadata_sample_file <- "./data_sample.metadata/metadata_samples_LLTK1toLLTK8_SalishSea.txt"
rubias_juvs_file <- "./analysis_rubias/analyses_manuscript_gutcontents/results_rubias_3repgroups_juvenilesalmon_byseason_posteriorprob.txt"
rubias_ads_file <- "./analysis_rubias/analyses_manuscript_gutcontents/results_rubias_3repgroups_adultsalmon_byseason_posteriorprob.txt"

##############################################################################
# Read in the data

meta_df <- read.delim(metadata_sample_file, sep = "\t", header = TRUE)
head(meta_df)

post_juv_df<- read.delim(rubias_juvs_file, sep = "\t", header = TRUE)
head(post_juv_df)

post_ads_df<- read.delim(rubias_ads_file, sep = "\t", header = TRUE)
head(post_ads_df)

#####################################################################################
# Process the data - juvenile salmon

# Combine the rubias results and sample metadata
juv_df <- dplyr :: left_join(post_juv_df, meta_df, by = c("indiv" = "genetics_id"))%>%
  filter(herring_length_mm != "NA") %>%
  mutate(herring_length_cm = herring_length_mm/10)

# Plot the juvenile salmon data
# Subset the data into layers for plotting
janfeb_juvs <- juv_df %>%
  filter(repunit == "Jan_Feb")

marapr_juvs<- juv_df %>%
  filter(repunit == "Mar_Apr")

plot1 <- ggplot( ) + 
  geom_jitter(data = marapr_juvs, aes(salmon_length_cm, herring_length_cm, shape = repunit, color = sum_post), width = 1, height = 1, size = 4, alpha = 0.7)+
  geom_jitter(data = janfeb_juvs, aes(salmon_length_cm, herring_length_cm, shape = repunit, color = sum_post), width = 1, height = 1, size = 4, alpha = 0.7)+
  facet_wrap(~mixture_collection)+
  xlim(0,25)+
  ylim(0,10)+
  theme_bw()+
  xlab("Salmon length (cm)")+
  ylab("Herring length (cm)")+
  scale_shape_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr"), values = c(16,17))+
  scale_color_viridis(name = "sum of posterior probability", begin = 1, end = 0, limits = c(0.6,1))

plot1


## Linear regression on the results

# For all juveniles
linearMod1 <- lm(herring_length_cm ~ salmon_length_cm, data= juv_df)  # build linear regression model on full data
print(linearMod1)
summary(linearMod1) #Adjusted R-squared:  0.009523 
#F-statistic: 1.692 on 1 and 71 DF,  p-value: 0.1975

# For juveniles by season

spring_juvs <- juv_df %>%
  filter(mixture_collection == "spring")

summer_juvs <- juv_df %>%
  filter(mixture_collection == "summer")

linearMod2 <- lm(herring_length_cm ~ salmon_length_cm, data= spring_juvs)  # build linear regression model on full data
print(linearMod2)
summary(linearMod2) #Adjusted R-squared:  -0.02317 
#F-statistic: 0.02607 on 1 and 42 DF,  p-value: 0.8725

linearMod3 <- lm(herring_length_cm ~ salmon_length_cm, data= summer_juvs)  # build linear regression model on full data
print(linearMod3)
summary(linearMod3) #Adjusted R-squared:  -0.03557 
#F-statistic: 0.03834 on 1 and 27 DF,  p-value: 0.8462 


#####################################################################################
# Process the data - adult salmon

# Combine the rubias results and sample metadata
adult_df <- dplyr :: left_join(post_ads_df, meta_df, by = c("indiv" = "genetics_id")) %>%
  filter(herring_length_mm != "NA") %>%
  mutate(herring_length_cm = herring_length_mm/10)


head(adult_df)

# make a plot

janfeb_ads <- adult_df %>%
  filter(repunit == "Jan_Feb")

marapr_ads<- adult_df %>%
  filter(repunit == "Mar_Apr")

may_ads<- adult_df %>%
  filter(repunit == "May_Jun")


plot2 <- ggplot( ) + 
  geom_jitter(data = marapr_ads, aes(salmon_length_cm, herring_length_cm, shape = repunit, color = sum_post), width = 1, height = 1, size = 4, alpha = 0.7)+
  geom_jitter(data = janfeb_ads, aes(salmon_length_cm, herring_length_cm, shape = repunit, color = sum_post), width = 1, height = 1, size = 4, alpha = 0.7)+
  geom_jitter(data = may_ads, aes(salmon_length_cm, herring_length_cm, shape = repunit, color = sum_post), width = 1, height = 1, size = 4, alpha = 0.7)+
  facet_wrap(~mixture_collection)+
  xlim(0,100)+
  ylim(0,25)+
  theme_bw()+
  xlab("Salmon length (cm)")+
  ylab("Herring length (cm)")+
  scale_shape_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = c(16,17,18))+
  scale_color_viridis(name = "sum of posterior probability", begin = 1, end = 0, limits = c(0.6,1))

plot2



## Linear regression on the results

# For all adults
linearMod4 <- lm(herring_length_cm ~ salmon_length_cm, data= adult_df)  # build linear regression model on full data
print(linearMod4)
summary(linearMod4) #Adjusted R-squared:  0.006495 
#F-statistic: 3.321 on 1 and 354 DF,  p-value: 0.06924


# For adults by season
spring_ads <- adult_df %>%
  filter(mixture_collection == "spring")

summer_ads <- adult_df %>%
  filter(mixture_collection == "summer")

winter_ads <- adult_df %>%
  filter(mixture_collection == "winter")

linearMod5 <- lm(herring_length_cm ~ salmon_length_cm, data= spring_ads)  # build linear regression model on full data
print(linearMod5)
summary(linearMod5) #Adjusted R-squared:  0.03731 
#F-statistic: 10.34 on 1 and 240 DF,  p-value: 0.00148

linearMod6 <- lm(herring_length_cm ~ salmon_length_cm, data= summer_ads)  # build linear regression model on full data
print(linearMod6)
summary(linearMod6) #Adjusted R-squared:  -0.01622 
#F-statistic: 0.01054 on 1 and 61 DF,  p-value: 0.9186

linearMod7 <- lm(herring_length_cm ~ salmon_length_cm, data= winter_ads)  # build linear regression model on full data
print(linearMod7)
summary(linearMod7) #Adjusted R-squared:  -0.01759 
#F-statistic: 0.1356 on 1 and 49 DF,  p-value: 0.7143


#####################################################################################
# make multiplot
final_plot <- plot_grid(
  plot1, plot2,
  labels = c("A", "B"),
  nrow = 2,
  align="v")

final_plot
#####################################################################################
# Save the plot to pdf file

ggsave("./analysis_fork.lengths/multiplot_compare_rubias_to_forklengths.pdf", final_plot)

