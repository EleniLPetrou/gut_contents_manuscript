# The purpose of this script is to conduct an AMOVA in R using a genepop file as input 

########################################################################################
# Load libraries
library(poppr)
library(pegas)
library(adegenet)
library(tidyverse)
library(genepopedit)

########################################################################################
# Specify the file names for the input files
genepop_file <- "./data_genepop.files/genepop_reference_7loci_20200710_filt.gen"
metadata_file <- "./data_sample.metadata/metadata_baseline.txt"

########################################################################################
# Read in genetic data and hierarchical strata files. 
data_all_loci <-read.genepop(genepop_file, ncode = 3)
levels(pop(data_all_loci))

# Change the format of the genepop file so that it is a dataframe that a human can work with
my_data_flat <- genepop_flatten(genepop_file)
my_data_flat[1:10, 1:10]

# Read in a dataframe containing information about each population's Julian date of sampling

julian_data <-read.delim(metadata_file)

julian_data <- julian_data %>%
  select(-code) #get rid of the code column that was used for making maps

########################################################################################
# Step #2: build a strata df that contains metadata for each individual

my_cols <- my_data_flat %>%
  select(SampleID, Population)

my_strata <- left_join(my_cols, julian_data, by = c("Population" = "population"))
head(my_strata)

########################################################################################
#  Step #3: assign those levels as strata to the genind object
strata(data_all_loci) <- my_strata

# Change the format of your data from a genind to a genclone object
my_data <- as.genclone(data_all_loci)
my_data

# View some hierarchical levels in the data set
table(strata(my_data, ~reporting_group))

########################################################################################
#  Step #4: Conduct the AMOVA
# Analysis of Molecular Variance (AMOVA) is a method of estimating population differentiation 
# directly from molecular data and testing hypotheses about such differentiation. 
# 1- Null hypothesis: samples are taken from a single panmictic population, or there is no difference in allele frequencies among the populations.

# Some relevant notes from the poppr package:
# "The implementation of AMOVA in poppr requires two very basic components: 
# (1) A distance matrix derived from the data and 
# (2) a separate table used to partition the data into different stratifications.
# The distance matrix can be calculated using any distance as long as it is euclidean.


# Because the molecular data consist of Euclidean distances derived from vectors of 1s and 0s, 
# the data are unlikely to follow a normal distribution. A null distribution is therefore 
# computed by resampling of the data (Excoffier, et al. 1992). In each permutation, each 
# individual is assigned to a randomly chosen population while holding the sample sizes constant. 
# These permutations are repeated many times, eventually building a null distribution. 
# Hypothesis testing is carried out relative to these resampling distributions.
# amova as implemented in pegas, performs a  hierarchical analysis of molecular variance as described in Excoffier et al. (1992). 
#This implementation accepts any number of hierarchical levels.

#The formula must be of the form d, ~ A/B/... where d is a genclone or genind object, and A, B, etc, 
# are the hierarchical levels from the highest to the lowest one. 
# Any number of levels is accepted, so specifying d ~ A will simply test for population differentiation.
# the poppr.amova function is a wrapper script for amova. It calculates a pairwise distance matrix by default."


# ADE4 AMOVAS
head(my_data)
amova_ade4 <- poppr.amova(my_data, ~reporting_group/location, within = FALSE, method = "ade4")
amova_ade4

amova.test <- randtest(amova_ade4, nrepet = 10000) # Test for significance
plot(amova.test)
amova.test
