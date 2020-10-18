# The purpose of this script is to conduct mixed stock analysis 
# using the R package rubias. 
# Data analyzed using "R version 3.6.1 (2019-07-05)"


###############################################################################################
# Load libraries
library(rubias)
library(tidyverse)


###############################################################################################
# Specify the input file names (containing reference samples in rubias format)
ref_file <- "./data_rubias.files/rubias_reference_7loci_20200710_filt_RubiasInput.txt"

# Specify the input file names (containing mixture samples in rubias format)
mix_file <- "./data_rubias.files/rubias_Marr15_Squa14_ElBy13_7loci_20200309_filt_RubiasInput.txt" 

###############################################################################################
# Specify the output file names:
output_mix_ests <- "./analysis_rubias/analyses_manuscript_baseline_validation/results_rubias_analysis_baseline_validation_3repgroups.txt"
output_post_probs <- "./analysis_rubias/analyses_manuscript_baseline_validation/results_rubias_baseline_validation_3repgroups_posteriorprob.txt"


###############################################################################################
# Read in data
ref_df <- read.delim(ref_file, header = TRUE, sep = " ", stringsAsFactors = FALSE, na.strings = "NA")
mix_df <- read.delim(mix_file, header = TRUE, sep = " ",  stringsAsFactors = FALSE, na.strings = "NA")

###############################################################################################
# Perform Genetic Mixture Analysis
#This is done with the infer_mixture function

mix_est <- infer_mixture(reference = ref_df, 
                         mixture = mix_df, 
                         gen_start_col = 5,
                         method = "MCMC",
                         reps = 10000, burn_in = 1000)

# The results comes are list of four tidy data frames that you can peek at:
lapply(mix_est, head)

# Explanation of terms:
# 1. mixing_proportions: the mixing proportions. The column pi holds the estimated mixing proportion for each collection.

#2. indiv_posteriors: this holds, for each individual, the posterior means of group membership 
# in each collection. Column PofZ holds those values. Column log_likelihood holds the log of the 
# probability of the individuals genotype given it is from the collection. Also included are n_non_miss_loci 
# and n_miss_loci which are the number of observed loci and the number of missing loci at the individual. 
# A list column missing_loci contains vectors with the indices (and the names) of the loci that are missing 
# in that individual. It also includes a column z_score which can be used to diagnose fish that don't belong 
# to any samples in the reference data base (see below).

#3. mix_prop_traces: MCMC traces of the mixing proportions for each collection. 
#You will use these if you want to make density estimates of the posterior distribution of 
#the mixing proportions or if you want to compute credible intervals.

###############################################################################################
# Aggregate collections into reporting units:

# Dataframe summarizing the estimated mixture proportions:
rep_mix_est <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # add mixing proportions over collections in the repunit

head(rep_mix_est)


# Dataframe summarizing the most likely reporting unit of origin for each individual, 
# by summing across the posterior probabilities across collections in a reporting group. 

options(scipen = 999) # turn off scientific notation

post_probs_df <- mix_est$indiv_posteriors %>%
  select(-missing_loci)  %>%
  group_by(indiv, repunit) %>%
  mutate(sum_post = sum(PofZ))%>%
  ungroup() %>%
  group_by(indiv) %>%
  top_n(1, PofZ) # this is just to get the top reporting group for each sample


###############################################################################################
# QC step : Assess whether individuals are not from any of the reference populations

# According to Eric Anderson via https://github.com/eriqande/rubias: 
# "Aberrantly low values of the genotype log-likelihood can indicate that there is something wrong. However, the raw likelihood that you get will depend 
# on the number of missing loci, etc. rubias deals with this by computing a z-score for each fish. 
# The Z-score is the Z statistic obtained from the fish's log-likelihood (by subtracting from it the 
# expected log-likelihood and dividing by the expected standard deviation). Here, we will look at the z-score computed for each fish to the population with the highest posterior. 
# (It is worth noting that you would never want to use the z-score to assign fish to different populations-it 
# is only there to decide whether it looks like it might not have actually come from the population that 
# it was assigned to, or any other population in the reference data set.)"


#If everything is ok, the z-scores we see will be roughly normally distributed.
normo <- tibble(z_score = rnorm(1e06))

# get the maximum-a-posteriori population for each individual
map_rows <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")



###############################################################################################
# Compute Credible Intervals for the mixture proportion estimates

# First, check how many MCMC sweeps were done:
nsweeps <- max(mix_est$mix_prop_traces$sweep)
nsweeps

# Discard the first 200 sweeps as burn-in
test <- mix_est$mix_prop_traces
head(test)

trace_subset <- mix_est$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(mixture_collection, repunit, sweep) %>%
  summarise(repprop = sum(pi))

head(trace_subset)

# now we can plot those traces and see what is up:
ggplot(trace_subset, aes(x = repprop, colour = repunit)) +
  geom_density()+
  facet_wrap(~mixture_collection)+
  theme_bw()

# next, let's take that information and calculate the 95% credible intervals from it
credible_interval_df <- trace_subset %>%
  group_by(mixture_collection, repunit) %>%
  summarise(loCI = round(quantile(repprop, probs = 0.05), 3),
            hiCI = round(quantile(repprop, probs = 0.95),3))

head(credible_interval_df)

## Next, add the credible interval df to the plotting_df, for later use
plotting_df <- left_join(rep_mix_est, credible_interval_df, by = c("mixture_collection", "repunit"))

plotting_df$repprop_round<- round(plotting_df$repprop, 3) #only look to 3 decimal places

head(plotting_df)


###############################################################################################
# Plot the results

my_cols <- c("#2b83ba","#abdda4",  "#fdae61") #some colors

mylabels = c(ElBy13 = "Elliot Bay, April 2013", Marr15=  "Cherry Point, May 2014", Squa14 = "Squaxin Pass, February 2014")

# order the label factor in a specific way

plotting_df$mixture_collection <- factor(plotting_df$mixture_collection, levels = c("Squa14", "ElBy13", "Marr15"))


# unstacked barplots, with credible intervals

assessment_plot <- ggplot(data = plotting_df, 
                          aes(x = "", 
                              y = repprop,
                              ymin= loCI, 
                              ymax= hiCI,
                              fill = repunit)) +
  geom_bar( stat = "identity", alpha = 0.8, position = "dodge")+
  geom_errorbar(position = position_dodge(), colour="dark grey") +
  ylab("Estimated mixing proportion") +
  xlab ("") +
  facet_wrap(~mixture_collection, labeller=labeller(mixture_collection = mylabels))+
  theme_classic()+
  scale_fill_manual(name = "Reporting group",labels = c("Jan-Feb", "Mar-Apr", "May"), values = my_cols)

assessment_plot 

###############################################################################################
# Save output of analyses to a text file

write.table(plotting_df, output_mix_ests, quote = FALSE, sep = "\t", row.names = FALSE)

write.table(post_probs_df, output_post_probs, quote = FALSE, sep = "\t", row.names = FALSE)


