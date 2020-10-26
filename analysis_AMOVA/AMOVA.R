# The purpose of this script is to conduct an AMOVA in R using a genotypic data in vcf format as input
# Load libraries
library(vcfR)
library(poppr)
library(pegas)
library(tidyverse)


######################################################################################
# Specify the names of data files used
vcf_fileName <- "data_vcf/batch_1_firstsnp_GCA900700415_mapq20_salish_7loci.recode.vcf"
metadata_fileName <- "data_sample.metadata/metadata_baseline.txt"

######################################################################################
# Specify the names of output file
out_fileName <- "analysis_AMOVA/amova_result_7loci.txt"

######################################################################################
# Read in a dataframe containing information about each population's Julian date of sampling
julian_data <-read.delim(metadata_fileName)
levels(julian_data$population)

# Read in the genotype data (with vcfR) and save it as a dataframe 
my_vcf <-  read.vcfR(vcf_fileName)

# Transform the genotype data into a genind object (intermediate step on the way to genclone object)
my_genind <- vcfR2genind(my_vcf)

remove(my_vcf)# you don't need the vcf in memory anymore


######################################################################################
# Do some data processing, in preparation for the analyses

# Save a vector of unique locus names
my_loci<- unique(my_genind$loc.fac)
length(my_loci) #check the number of loci

# Create a vector of individual names
name_vec <- indNames(my_genind)
head(name_vec)

# Create a vector of population names
pop_vec <- sapply(strsplit(name_vec, "_"), `[`, 1)
head(pop_vec)
# Assign those population names to your genind object
pop(my_genind) <- pop_vec
my_genind@pop #check

# bind the vectors together into a df
my_df <- as.data.frame(cbind(name_vec, pop_vec))
head(my_df)

######################################################################################
# Step #2: build a strata df that contains metadata for each individual

my_strata <- left_join(my_df, julian_data, by = c("pop_vec" = "population"))
head(my_strata)

######################################################################################
# Step #3: assign the strata df as strata to the genind object

strata(my_genind) <- my_strata
strata(my_genind)

# View some hierarchical levels in the data set
table(strata(my_genind, ~reporting_group))

######################################################################################
# ADE4 AMOVAS #
# AMOVA is a classical method of assessing population differentiation by evaluating where 
# the most variation exists in a hierarchical population structure (Excoffier et al., 1992).
# The ADE4 Amova

sink(out_fileName) #start writing output to file

amova_ade4 <- poppr.amova(my_genind, ~reporting_group/location, within = FALSE, method = "ade4")
amova_ade4

amova.test <- randtest(amova_ade4, nrepet = 999) # Test for significance
plot(amova.test)
amova.test

sink() # Stop writing to file


###pegas AMOVAS#######################################################
# The pegas Amova report phi statistics in an easy to understand way. 
# I like runing both AMOVAs and cross-validating them. 
# To run an AMOVA in pegas, you need two things:
  
#1. a distance matrix among samples
#2.  a data frame with hierarchical population strata
herring_dist  <- dist(my_genind) #  Calculate euclidean distance btwn pops
herring_stra  <- strata(my_genind)

# ready for the AMOVA!
amova_pegas <- pegas::amova(herring_dist ~ reporting_group/location, data = herring_stra, nperm = 999)
amova_pegas 
