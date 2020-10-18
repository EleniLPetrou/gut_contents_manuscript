# The purpose of this script is to take a file of genotypes saved in genepop format
# and convert it into rubias format. 
# Script by Eleni Petrou 20201016, using "R version 3.6.1 (2019-07-05)"


# To do this, you need to install "genepopedit"
devtools::install_github("rystanley/genepopedit")

############################################################################
# Load the necessary libraries
library (tidyverse)
library(genepopedit)

############################################################################
# Part 1: Create a rubias mixture file for the gut content samples

# Filename of mixed samples in genepop format
mix_genepop1 <- "./data_genepop.files/genepop_gutsamples_7loci_20200710_filt.gen"

# Set output path
outpath1 <- "./data_rubias.files/rubias_gutsamples_7loci_20200710_filt_"

# Count the number of individuals in your data set
n_mix1 <-length(genepop_detective(mix_genepop1, variable="Inds"))
n_mix1

# Populate the repunit column of rubias file
mix_repunit1 <- rep("NA", n_mix1)

# Make the rubias file and save it to disk
# sampletype - specifies whether your rubias input file is a 
#reference (baseline) or mixture file. Can only be one of "reference" or "mixture"

genepop_rubias(genepop = mix_genepop1,  
               sampletype = "mixture", 
               repunits = mix_repunit1, 
               path = outpath1 )


############################################################################
# Part 2: Create a rubias mixture file for the GSI validation samples 
# (independent samples collected at spawning grounds and not used for locus discovery)

# Filename of samples in genepop format
mix_genepop2 <- "./data_genepop.files/genepop_Marr15_Squa14_ElBy13_7loci_20200309_filt.gen"

# Set output path
outpath2 <- "./data_rubias.files/rubias_Marr15_Squa14_ElBy13_7loci_20200309_filt_"

# Count the number of individuals in your data set
n_mix2 <-length(genepop_detective(mix_genepop2, variable="Inds"))
n_mix2

# Populate the repunit column of rubias file
mix_repunit2 <- rep("NA", n_mix2)

# # Make the rubias file and save it to disk
# sampletype - specifies whether your rubias input file is a 
#reference (baseline) or mixture file. Can only be one of "reference" or "mixture"

genepop_rubias(genepop = mix_genepop2,  
               sampletype = "mixture", 
               repunits = mix_repunit2, 
               path = outpath2 )


############################################################################
# Part 3: Create a rubias reference file using the genetic baseline samples

# File name of genepop file containing the baseline samples
ref_genepop <- "./data_genepop.files/genepop_reference_7loci_20200710_filt.gen"

# Set output path
outpath3 <- "./data_rubias.files/rubias_reference_7loci_20200710_filt_"

#File name of txt file containing repunit information for each baseline sample
ref_metadata <- "./data_rubias.files/rubias_baseline_repunits.txt"
repunit_df <- read.delim(ref_metadata) #read in the metadata

# Count the number of individuals in your data set
n_ref <-length(genepop_detective(ref_genepop, variable="Inds"))
n_ref

# Populate the repunit column of rubias file
ref_repunit <- repunit_df$repunit3

# # Make the rubias file and save it to disk
# sampletype - specifies whether your rubias input file is a 
#reference (baseline) or mixture file. Can only be one of "reference" or "mixture"

genepop_rubias(genepop = ref_genepop,  
               sampletype = "reference", 
               repunits = ref_repunit, 
               path = outpath3 )


