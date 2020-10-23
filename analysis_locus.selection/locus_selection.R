# locus selection script
# The purpose of this script is to identify loci that are divergent between herring populations
# "R version 3.6.1 (2019-07-05)"


# Load libraries
library(adegenet)
library(hierfstat)
library(tidyverse)
library(viridis)
library(vcfR)


######################################################################################
# input fileName (vcf)
fileName<- "./data_vcf/batch_1_firstsnp_GCA900700415_mapq20_salish.recode.vcf"

######################################################################################
#output filename (txt)
output_fileName<- "./analysis_locus.selection/taqman_candidates.txt"
######################################################################################
# Read in the data (with vcfR)
my_vcf <-  read.vcfR(fileName)

# Transform the data into a genind object, that hierfstat and adegenet can use
my_genind <- vcfR2genind(my_vcf)

######################################################################################
# Do some data processing, in preparation for the analyses
# Save a vector of unique locus names

my_loci<- unique(my_genind$loc.fac)
length(my_loci) #check the number of loci

# Create a vector of individual names
name_vec <- indNames(my_genind)
head(name_vec) #check

# Create a vector of population names
pop_vec <- sapply(strsplit(name_vec, "_"), `[`, 1)
head(pop_vec) #check

# Assign those population names to your genind object
pop(my_genind) <- pop_vec
my_genind@pop #check
indNames(my_genind)

##################################################################################################
#  Part 1: Conduct DAPC

# Find clusters of populations in the data.
# Then specify the maximum number of PCAs to retain (by looking at the graph). There is no reason to discard PCAs
# other than computational time, so keep all of them.
# Finally, look at the graph output: Value of BIC versus number of clusters; the optimal number of clusters
# is the minimum BIC score (at the elbow of the graph).
# DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components.
popn <- length(levels(my_genind@pop)) 

dapc_all <- dapc(my_genind,my_genind$pop,n.pca=350,n.da=(popn-1)) ##Retain all, then identify optimal number by optim.a.score

# Look at the optim_a score and choose the optimal number of PCS based on the highest value
test_a_score <- optim.a.score(dapc_all)
test_a_score$best #best score

dapc_all <- dapc(my_genind,my_genind$pop,n.pca=test_a_score$best,n.da=popn) ##Here, 64 PC's is the optimal number
dapc_all

# Save the output of the individual coordinates from the DAPC as a dataframe for the first six discriminant axes. 
dapc_df <- as.data.frame((dapc_all$ind.coord[ , 1:6]))
head(dapc_df)

# Add the sample id and population id metadata to this dataframe
dapc_df$Sample <- name_vec
dapc_df$Population <- pop_vec

# Plot the DAPC results
ggplot(data = dapc_df, aes(x= LD1, y= LD2)) + 
  geom_point(aes(color = Population))+
  theme_bw()

#######################################################################################
# Assess which alleles pull apart the DAPC clusters using the command loadingplot.
# Variable contributions are stored in the var.contr slot of a dapc object
#Along Axis 1
loadingplot(dapc_all$var.contr, axis = 1, lab.jitter = 1)

# Along Axis 2
loadingplot(dapc_all$var.contr, axis = 2, lab.jitter = 1)

# Save the DAPC loadings of the first three discriminant axes to a dataframe
loadings_df <- as.data.frame(dapc_all$var.contr) %>%
  rownames_to_column("loc_all") %>%
  separate(loc_all, c("ID", "Allele"), "\\.") %>% #split this idiotic string into two useful variables, the locus name and the allele
  filter(Allele == 0) %>% # save only one of the alleles to avoid redundancy
  select(ID, LD1, LD2, LD3) %>% #output specific columns
  arrange(desc(LD1)) #sort the dataframe so the locus with the largest value along DA1 is at the top

head(loadings_df)

# Add some important locus metadata to the DAPC loadings.
locus_metadata <- as.data.frame(my_vcf@fix) %>%
  select(CHROM, POS, ID, REF, ALT)

head(locus_metadata)

tidy_df<- left_join(loadings_df, locus_metadata, by = "ID") 
head(tidy_df)

# From each chromosome, select the top 20 loci that have the highest loading scores

LD1_tops <- tidy_df %>%
  group_by(CHROM) %>%
  top_n(n = 20, wt= LD1)


##################################################################################################
#  Part 2: Calculate global FST per locus
# Calculate a suite of basic population genetics metrics in your dataset using the hierfstat package

my_stats <- basic.stats(my_genind)

# Get Fst per locus and save the output to a df

stats_df <- my_stats$perloc %>%
  rownames_to_column("ID")

head(stats_df)

# join output with tidy_df

final_df<- left_join(tidy_df, stats_df, by = "ID") 
head(final_df)

# Retain loci in final_df if they had a global FST > 0.1, Ho>0.1, and were ramked as the top 20 loci in their chromosome for LD1 or LD2

taqmen_set1 <- final_df %>%
  filter(Fst > 0.10, Ho > 0.1)%>%
  group_by(CHROM) %>%
  top_n(n = 20, wt = LD1)

taqmen_set2 <- final_df %>%
  filter(Fst > 0.10, Ho > 0.1)%>%
  group_by(CHROM) %>%
  top_n(n = 20, wt = LD2)
  
#######################################################################################
# Save the data to text files
write.table(taqmen_set1, file = output_fileName, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE)
