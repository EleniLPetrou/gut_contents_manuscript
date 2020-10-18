# This script was written by Eleni Petrou on 20201015 using "R version 3.6.1 (2019-07-05)"
# Its purpose is to conduct a PCA using a genepop file and sample metadata files as input data

############################################################################################### 
# Load the necessary libraries
library(adegenet)
library (tidyverse)
library(cowplot)
library(viridis)
library(scales)

############################################################################################### 
# Specify the file names
input_genepop <- "./data_genepop.files/genepop_refandmixcombo_7loci_20200710_filt.gen"
metadata_baseline_file <- "./data_sample.metadata/metadata_baseline.txt"
metadata_sample_file <- "./data_sample.metadata/metadata_samples_LLTK1toLLTK8_SalishSea.txt"

############################################################################################### 
# Read in the files

# First, read in your data as a genepop file. The file can be delimited by tabs or spaces but there must abe a comma after each individual. 
# Specify how many characters code each allele with ncode. 
geno_genind <-read.genepop(input_genepop, ncode = 3)
(summary(geno_genind))

# Read in the sample metadata for plotting
metadata_baseline_df <- read.delim(metadata_baseline_file, sep = "\t", header = TRUE)

metadata_sample_df <- read.delim(metadata_sample_file, sep = "\t", header = TRUE, na.strings = "NA")

# Check to make sure that fish lengths are stored as numeric data
class(metadata_sample_df$salmon_length_cm)
class(metadata_sample_df$herring_length_mm)

############################################################################################### 
# Conduct PCA

# Replace missing data by the mean allele frequencies
geno_genind_scaled <- scaleGen(geno_genind, NA.method="mean")

# PCA
pca_A <- dudi.pca(geno_genind_scaled,cent=TRUE,scale=FALSE,scannf=FALSE,nf=3)
summary(pca_A) # The Projected inertia columns store information on the percent variance explained by each axis

# Save the output of the PCA to a dataframe for plotting and 
# Reorganize the dataframe such that you get sample and population names
A_df <- pca_A$li %>% 
  tibble :: rownames_to_column(var = "sample") %>%
  dplyr :: mutate(sample2 = sample) %>%
  tidyr :: separate(sample2, c("population", "temp"), "_") %>%
  dplyr :: select(-temp)
  

# combine the PCA results df with the sample metadata to make some nice plots
# For baseline metadata
head(A_df2)
head(metadata_baseline_df)

PCA_df <- dplyr :: left_join(A_df, metadata_baseline_df, by = "population")
head(PCA_df)

############################################################################################### 
# Prepare the data for plotting

# First, filter the data, so you can add it as individual layers
# subset the PCA_df so you can plot the baseline samples separately from the gut-content samples
baseline_samples <- PCA_df %>%
  dplyr :: filter(population != "will" , population != "josh" ) %>%
  select(-code, -dummy_date)

gut_samples <- PCA_df %>%
  dplyr :: filter( population == "josh" | population == "will" ) %>% # samples collected by josh and will are the gut contents samples
  dplyr :: select(sample, Axis1, Axis2, Axis3, population) %>%
  dplyr :: left_join(metadata_sample_df, by = c( "sample" = "genetics_id"))%>% #add metadata to each sample
  mutate(salmon_class = if_else(salmon_length_cm > 25, "adult", "juvenile")) #classify each sample as being from a juvenile or adult salmon


########################################################################
# Plot the PCA

# set the breaks and labels for your color scheme
mybreaks=c(0, 30, 60, 90, 120, 150)
mylabels = c("January", "February", "March", "April", "May", "June")


plotA <- ggplot() +
  geom_point(data = baseline_samples, aes(x = Axis1, y = Axis2, color = julian_date),  size = 2, alpha = 0.6, shape = 16)+
  geom_point(data = gut_samples, aes(x = Axis1, y = Axis2, shape = population), size = 4, alpha = 0.6)+
  scale_color_viridis(option="plasma", name="Collection date of \nspawning herring", labels = mylabels,breaks = mybreaks, begin = 0, end = 0.9)+
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
 

plotA

########################################################################
# Write out the PCA dataframes to tab-delimited text files

write.table(baseline_samples, 
            "./analysis_PCA/results_PCA_baseline_samples.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(gut_samples, 
            "./analysis_PCA/results_PCA_gut_samples.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

