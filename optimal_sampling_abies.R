# A GenAlEx formatted excel sheet (.xlsx) is the required input file.

# Set working directory, where the input file is. On RStudio you can achieve 
#   this by navigating to "Session" -> "Set Working Directory" -> 
#   "Choose Directory"

# Warning of entirely non-type individuals when importing datasets

# After running some of the commands the following warnings appear: 
# "In validityMethod(object) :
# @tab does not contain integers; as of adegenet_2.0-0, numeric values are no longer used"
# They are triggered because of replacing missing data with the mean value. 
# In this case they are harmless.
# https://groups.google.com/forum/#!topic/poppr/F-HImtnMrA8

# Due to a bug, package hierfstat is incapable of calculating allelic richness for 
# single-locus datasets.
# A bug report has been sent: https://github.com/jgx65/hierfstat/issues/25
# The workaround implemented in this script, is to import a new dataset (so a new excel 
# sheet has to be created) where the single - most polymorphic locus is duplicated. 
# In that way, mean measures of Ar are unchanged.
# I have created such a sheet for locus "SFb4", which according to the 
# common LifeGenMon dataset is the most polymorphic locus for Abies. 

# Install (if needed) & load required libraries
if (!require("adegenet")) install.packages("adegenet")
library(adegenet)
if (!require("popprxl")) install.packages("popprxl")
library(popprxl)
if (!require("hierfstat")) install.packages("hierfstat")
library(hierfstat)
if (!require("reshape2")) install.packages("reshape2")
library(reshape2)
if (!require("ggplot2")) install.packages("ggplot2") 
library(ggplot2)
if (!require("dplyr")) install.packages("dplyr") 
library(dplyr)

# Load dataset
# A GenAlEx formatted excel sheet is the required input
obj <- read.genalexcel(
  "LGM_DE_SI_GR_final.xlsx",   # name of excel file
  sheet = "Abies",             # name of sheet where the genotypes reside
  genclone = F) 

#nSSRs <- c("SF1", "NFF3", "Aag01", "NFH15", "NFF7", "SFb4")
#EST_SSRs <- c("Aat06", "Aat11", "Aat15", "Aat01", "Aat04")
#obj <- obj[loc = EST_SSRs]

loci <- sort(nAll(obj)) # vector containing loci from least to
# most polymorphic according to the pooled dataset from all countries and pops
loci <- names(loci) # transform named vector to vector of names
most_poly_locus <- loci[length(loci)]

# There is a bug in the current versions of hierfstat (on cran & development version)
# which prevents the package from calculating Ar for single-locus datasets.
# The workaround implemented here, is to import a new dataset where the single locus
# is duplicated, so mean measures of Ar are unaffected.

paste("A sheet in the excel input file needs to be created for locus", most_poly_locus, sep = " ")

# Import single-locus dataset
obj_fix <- read.genalexcel(
  "LGM_DE_SI_GR_final.xlsx",     # name of excel file
  sheet = "Abies-SFb4",     # name of sheet of single-locus dataset
  genclone = F)

# Real-world datasets are expected to contain missing data which might introduce bias in 
# simulations. With the following function, missing data are replaced with mean values.
# Comment next line to leave missing data as they are.
obj <- missingno(obj, type = "mean")

pop <- "DE_Regen" # select pop to analyze

replic_num <- 100   # set number of replications

# Simulations ####
system.time({
  
  # Seperate dataset by pop & select pop to analyze
  obj_list <- seppop(obj) # separate pops
  obj <- obj_list[[pop]] 
  
  # Set sample size
  if(pop == "DE_Adult"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250")
  }else if (pop == "GR_Adult"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250")
  }else if (pop == "SL_Adult"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","249")
  }else if (pop == "DE_Regen"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200")
  }else if (pop == "GR_Regen"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200")
  }else if (pop == "SL_Regen"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200")
  }else if (pop == "DE_Seed"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250","275", "300", "325", "350", "375", "400")
  }else if (pop == "GR_Seed"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250","275", "300", "325", "350", "375", "382")
  }else if (pop == "SL_Seed"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250","275", "300", "325", "350", "375", "400")
  }else{
    print("Unknown population")
  }
  
  data <- list()
  for(i in 1:length(loci)){
    data[[length(loci)+1-i]] <- obj[, loc = loci[i:length(loci)]]
  }
  
  # Save the highest sample size 
  high_samp_size <- samp_size[length(samp_size)]# added because of hierfstat bug 
  
  # sim_data_17 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_17 <- list()
  for(i in samp_size){ 
    sim_data_17[[i]] <-
      replicate (replic_num, data[[17]][sample(1:nrow(data[[17]]$tab), 
                                               i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_17[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_17 <- melt(Hobs_means)
  colnames(Hobs_means_df_17) <- c("value", "replic", "samp_size")
  Hobs_means_df_17$samp_size <- 
    factor(Hobs_means_df_17$samp_size, levels = unique(
      as.character(Hobs_means_df_17$samp_size)))
  
  Hobs_means_df_17$marker_num <- "17 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_17 <- melt(Hexp_means)
  colnames(Hexp_means_df_17) <- c("value", "replic", "samp_size")
  Hexp_means_df_17$samp_size <- 
    factor(Hexp_means_df_17$samp_size, levels = unique(
      as.character(Hexp_means_df_17$samp_size)))
  
  Hexp_means_df_17$marker_num <- "17 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_17[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_17 <- melt(ar_means)
  colnames(ar_means_df_17) <- c("value", "samp_size")
  ar_means_df_17$samp_size <- 
    factor(ar_means_df_17$samp_size, levels = unique(
      as.character(ar_means_df_17$samp_size)))
  
  ar_means_df_17$marker_num <- "17 markers"
  # sim_data_16 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_16 <- list()
  for(i in samp_size){ 
    sim_data_16[[i]] <-
      replicate (replic_num, data[[16]][sample(1:nrow(data[[16]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_16[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_16 <- melt(Hobs_means)
  colnames(Hobs_means_df_16) <- c("value", "replic", "samp_size")
  Hobs_means_df_16$samp_size <- 
    factor(Hobs_means_df_16$samp_size, levels = unique(
      as.character(Hobs_means_df_16$samp_size)))
  
  Hobs_means_df_16$marker_num <- "16 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_16 <- melt(Hexp_means)
  colnames(Hexp_means_df_16) <- c("value", "replic", "samp_size")
  Hexp_means_df_16$samp_size <- 
    factor(Hexp_means_df_16$samp_size, levels = unique(
      as.character(Hexp_means_df_16$samp_size)))
  
  Hexp_means_df_16$marker_num <- "16 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_16[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_16 <- melt(ar_means)
  colnames(ar_means_df_16) <- c("value", "samp_size")
  ar_means_df_16$samp_size <- 
    factor(ar_means_df_16$samp_size, levels = unique(
      as.character(ar_means_df_16$samp_size)))
  
  ar_means_df_16$marker_num <- "16 markers"
  # sim_data_15 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_15 <- list()
  for(i in samp_size){ 
    sim_data_15[[i]] <-
      replicate (replic_num, data[[15]][sample(1:nrow(data[[15]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_15[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_15 <- melt(Hobs_means)
  colnames(Hobs_means_df_15) <- c("value", "replic", "samp_size")
  Hobs_means_df_15$samp_size <- 
    factor(Hobs_means_df_15$samp_size, levels = unique(
      as.character(Hobs_means_df_15$samp_size)))
  
  Hobs_means_df_15$marker_num <- "15 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_15 <- melt(Hexp_means)
  colnames(Hexp_means_df_15) <- c("value", "replic", "samp_size")
  Hexp_means_df_15$samp_size <- 
    factor(Hexp_means_df_15$samp_size, levels = unique(
      as.character(Hexp_means_df_15$samp_size)))
  
  Hexp_means_df_15$marker_num <- "15 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_15[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_15 <- melt(ar_means)
  colnames(ar_means_df_15) <- c("value", "samp_size")
  ar_means_df_15$samp_size <- 
    factor(ar_means_df_15$samp_size, levels = unique(
      as.character(ar_means_df_15$samp_size)))
  
  ar_means_df_15$marker_num <- "15 markers"
  # sim_data_14 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_14 <- list()
  for(i in samp_size){ 
    sim_data_14[[i]] <-
      replicate (replic_num, data[[14]][sample(1:nrow(data[[14]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_14[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_14 <- melt(Hobs_means)
  colnames(Hobs_means_df_14) <- c("value", "replic", "samp_size")
  Hobs_means_df_14$samp_size <- 
    factor(Hobs_means_df_14$samp_size, levels = unique(
      as.character(Hobs_means_df_14$samp_size)))
  
  Hobs_means_df_14$marker_num <- "14 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_14 <- melt(Hexp_means)
  colnames(Hexp_means_df_14) <- c("value", "replic", "samp_size")
  Hexp_means_df_14$samp_size <- 
    factor(Hexp_means_df_14$samp_size, levels = unique(
      as.character(Hexp_means_df_14$samp_size)))
  
  Hexp_means_df_14$marker_num <- "14 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_14[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_14 <- melt(ar_means)
  colnames(ar_means_df_14) <- c("value", "samp_size")
  ar_means_df_14$samp_size <- 
    factor(ar_means_df_14$samp_size, levels = unique(
      as.character(ar_means_df_14$samp_size)))
  
  ar_means_df_14$marker_num <- "14 markers"
  # sim_data_13 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_13 <- list()
  for(i in samp_size){ 
    sim_data_13[[i]] <-
      replicate (replic_num, data[[13]][sample(1:nrow(data[[13]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_13[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_13 <- melt(Hobs_means)
  colnames(Hobs_means_df_13) <- c("value", "replic", "samp_size")
  Hobs_means_df_13$samp_size <- 
    factor(Hobs_means_df_13$samp_size, levels = unique(
      as.character(Hobs_means_df_13$samp_size)))
  
  Hobs_means_df_13$marker_num <- "13 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_13 <- melt(Hexp_means)
  colnames(Hexp_means_df_13) <- c("value", "replic", "samp_size")
  Hexp_means_df_13$samp_size <- 
    factor(Hexp_means_df_13$samp_size, levels = unique(
      as.character(Hexp_means_df_13$samp_size)))
  
  Hexp_means_df_13$marker_num <- "13 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_13[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_13 <- melt(ar_means)
  colnames(ar_means_df_13) <- c("value", "samp_size")
  ar_means_df_13$samp_size <- 
    factor(ar_means_df_13$samp_size, levels = unique(
      as.character(ar_means_df_13$samp_size)))
  
  ar_means_df_13$marker_num <- "13 markers"
  # sim_data_12 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_12 <- list()
  for(i in samp_size){ 
    sim_data_12[[i]] <-
      replicate (replic_num, data[[12]][sample(1:nrow(data[[12]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_12[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_12 <- melt(Hobs_means)
  colnames(Hobs_means_df_12) <- c("value", "replic", "samp_size")
  Hobs_means_df_12$samp_size <- 
    factor(Hobs_means_df_12$samp_size, levels = unique(
      as.character(Hobs_means_df_12$samp_size)))
  
  Hobs_means_df_12$marker_num <- "12 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_12 <- melt(Hexp_means)
  colnames(Hexp_means_df_12) <- c("value", "replic", "samp_size")
  Hexp_means_df_12$samp_size <- 
    factor(Hexp_means_df_12$samp_size, levels = unique(
      as.character(Hexp_means_df_12$samp_size)))
  
  Hexp_means_df_12$marker_num <- "12 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_12[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_12 <- melt(ar_means)
  colnames(ar_means_df_12) <- c("value", "samp_size")
  ar_means_df_12$samp_size <- 
    factor(ar_means_df_12$samp_size, levels = unique(
      as.character(ar_means_df_12$samp_size)))
  
  ar_means_df_12$marker_num <- "12 markers"
  # sim_data_11 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_11 <- list()
  for(i in samp_size){ 
    sim_data_11[[i]] <-
      replicate (replic_num, data[[11]][sample(1:nrow(data[[11]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_11[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_11 <- melt(Hobs_means)
  colnames(Hobs_means_df_11) <- c("value", "replic", "samp_size")
  Hobs_means_df_11$samp_size <- 
    factor(Hobs_means_df_11$samp_size, levels = unique(
      as.character(Hobs_means_df_11$samp_size)))
  
  Hobs_means_df_11$marker_num <- "11 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_11 <- melt(Hexp_means)
  colnames(Hexp_means_df_11) <- c("value", "replic", "samp_size")
  Hexp_means_df_11$samp_size <- 
    factor(Hexp_means_df_11$samp_size, levels = unique(
      as.character(Hexp_means_df_11$samp_size)))
  
  Hexp_means_df_11$marker_num <- "11 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_11[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_11 <- melt(ar_means)
  colnames(ar_means_df_11) <- c("value", "samp_size")
  ar_means_df_11$samp_size <- 
    factor(ar_means_df_11$samp_size, levels = unique(
      as.character(ar_means_df_11$samp_size)))
  
  ar_means_df_11$marker_num <- "11 markers"
  # sim_data_10 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_10 <- list()
  for(i in samp_size){ 
    sim_data_10[[i]] <-
      replicate (replic_num, data[[10]][sample(1:nrow(data[[10]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_10[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_10 <- melt(Hobs_means)
  colnames(Hobs_means_df_10) <- c("value", "replic", "samp_size")
  Hobs_means_df_10$samp_size <- 
    factor(Hobs_means_df_10$samp_size, levels = unique(
      as.character(Hobs_means_df_10$samp_size)))
  
  Hobs_means_df_10$marker_num <- "10 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_10 <- melt(Hexp_means)
  colnames(Hexp_means_df_10) <- c("value", "replic", "samp_size")
  Hexp_means_df_10$samp_size <- 
    factor(Hexp_means_df_10$samp_size, levels = unique(
      as.character(Hexp_means_df_10$samp_size)))
  
  Hexp_means_df_10$marker_num <- "10 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_10[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_10 <- melt(ar_means)
  colnames(ar_means_df_10) <- c("value", "samp_size")
  ar_means_df_10$samp_size <- 
    factor(ar_means_df_10$samp_size, levels = unique(
      as.character(ar_means_df_10$samp_size)))
  
  ar_means_df_10$marker_num <- "10 markers"  
  # sim_data_09 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_09 <- list()
  for(i in samp_size){ 
    sim_data_09[[i]] <-
      replicate (replic_num, data[[9]][sample(1:nrow(data[[9]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_09[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_09 <- melt(Hobs_means)
  colnames(Hobs_means_df_09) <- c("value", "replic", "samp_size")
  Hobs_means_df_09$samp_size <- 
    factor(Hobs_means_df_09$samp_size, levels = unique(
      as.character(Hobs_means_df_09$samp_size)))
  
  Hobs_means_df_09$marker_num <- "9 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_09 <- melt(Hexp_means)
  colnames(Hexp_means_df_09) <- c("value", "replic", "samp_size")
  Hexp_means_df_09$samp_size <- 
    factor(Hexp_means_df_09$samp_size, levels = unique(
      as.character(Hexp_means_df_09$samp_size)))
  
  Hexp_means_df_09$marker_num <- "9 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_09[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_09 <- melt(ar_means)
  colnames(ar_means_df_09) <- c("value", "samp_size")
  ar_means_df_09$samp_size <- 
    factor(ar_means_df_09$samp_size, levels = unique(
      as.character(ar_means_df_09$samp_size)))
  
  ar_means_df_09$marker_num <- "9 markers"
  # sim_data_08 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_08 <- list()
  for(i in samp_size){ 
    sim_data_08[[i]] <-
      replicate (replic_num, data[[8]][sample(1:nrow(data[[8]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_08[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_08 <- melt(Hobs_means)
  colnames(Hobs_means_df_08) <- c("value", "replic", "samp_size")
  Hobs_means_df_08$samp_size <- 
    factor(Hobs_means_df_08$samp_size, levels = unique(
      as.character(Hobs_means_df_08$samp_size)))
  
  Hobs_means_df_08$marker_num <- "8 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_08 <- melt(Hexp_means)
  colnames(Hexp_means_df_08) <- c("value", "replic", "samp_size")
  Hexp_means_df_08$samp_size <- 
    factor(Hexp_means_df_08$samp_size, levels = unique(
      as.character(Hexp_means_df_08$samp_size)))
  
  Hexp_means_df_08$marker_num <- "8 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_08[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_08 <- melt(ar_means)
  colnames(ar_means_df_08) <- c("value", "samp_size")
  ar_means_df_08$samp_size <- 
    factor(ar_means_df_08$samp_size, levels = unique(
      as.character(ar_means_df_08$samp_size)))
  
  ar_means_df_08$marker_num <- "8 markers"
  # sim_data_07 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_07 <- list()
  for(i in samp_size){ 
    sim_data_07[[i]] <-
      replicate (replic_num, data[[7]][sample(1:nrow(data[[7]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_07[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_07 <- melt(Hobs_means)
  colnames(Hobs_means_df_07) <- c("value", "replic", "samp_size")
  Hobs_means_df_07$samp_size <- 
    factor(Hobs_means_df_07$samp_size, levels = unique(
      as.character(Hobs_means_df_07$samp_size)))
  
  Hobs_means_df_07$marker_num <- "7 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_07 <- melt(Hexp_means)
  colnames(Hexp_means_df_07) <- c("value", "replic", "samp_size")
  Hexp_means_df_07$samp_size <- 
    factor(Hexp_means_df_07$samp_size, levels = unique(
      as.character(Hexp_means_df_07$samp_size)))
  
  Hexp_means_df_07$marker_num <- "7 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_07[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_07 <- melt(ar_means)
  colnames(ar_means_df_07) <- c("value", "samp_size")
  ar_means_df_07$samp_size <- 
    factor(ar_means_df_07$samp_size, levels = unique(
      as.character(ar_means_df_07$samp_size)))
  
  ar_means_df_07$marker_num <- "7 markers"
  # sim_data_06 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_06 <- list()
  for(i in samp_size){ 
    sim_data_06[[i]] <-
      replicate (replic_num, data[[6]][sample(1:nrow(data[[6]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_06[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_06 <- melt(Hobs_means)
  colnames(Hobs_means_df_06) <- c("value", "replic", "samp_size")
  Hobs_means_df_06$samp_size <- 
    factor(Hobs_means_df_06$samp_size, levels = unique(
      as.character(Hobs_means_df_06$samp_size)))
  
  Hobs_means_df_06$marker_num <- "6 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_06 <- melt(Hexp_means)
  colnames(Hexp_means_df_06) <- c("value", "replic", "samp_size")
  Hexp_means_df_06$samp_size <- 
    factor(Hexp_means_df_06$samp_size, levels = unique(
      as.character(Hexp_means_df_06$samp_size)))
  
  Hexp_means_df_06$marker_num <- "6 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_06[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_06 <- melt(ar_means)
  colnames(ar_means_df_06) <- c("value", "samp_size")
  ar_means_df_06$samp_size <- 
    factor(ar_means_df_06$samp_size, levels = unique(
      as.character(ar_means_df_06$samp_size)))
  
  ar_means_df_06$marker_num <- "6 markers"
  # sim_data_05 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_05 <- list()
  for(i in samp_size){ 
    sim_data_05[[i]] <-
      replicate (replic_num, data[[5]][sample(1:nrow(data[[5]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_05[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_05 <- melt(Hobs_means)
  colnames(Hobs_means_df_05) <- c("value", "replic", "samp_size")
  Hobs_means_df_05$samp_size <- 
    factor(Hobs_means_df_05$samp_size, levels = unique(
      as.character(Hobs_means_df_05$samp_size)))
  
  Hobs_means_df_05$marker_num <- "5 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_05 <- melt(Hexp_means)
  colnames(Hexp_means_df_05) <- c("value", "replic", "samp_size")
  Hexp_means_df_05$samp_size <- 
    factor(Hexp_means_df_05$samp_size, levels = unique(
      as.character(Hexp_means_df_05$samp_size)))
  
  Hexp_means_df_05$marker_num <- "5 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_05[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_05 <- melt(ar_means)
  colnames(ar_means_df_05) <- c("value", "samp_size")
  ar_means_df_05$samp_size <- 
    factor(ar_means_df_05$samp_size, levels = unique(
      as.character(ar_means_df_05$samp_size)))
  
  ar_means_df_05$marker_num <- "5 markers"
  
  # sim_data_04 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_04 <- list()
  for(i in samp_size){ 
    sim_data_04[[i]] <-
      replicate (replic_num, data[[4]][sample(1:nrow(data[[4]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_04[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_04 <- melt(Hobs_means)
  colnames(Hobs_means_df_04) <- c("value", "replic", "samp_size")
  Hobs_means_df_04$samp_size <- 
    factor(Hobs_means_df_04$samp_size, levels = unique(
      as.character(Hobs_means_df_04$samp_size)))
  
  Hobs_means_df_04$marker_num <- "4 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_04 <- melt(Hexp_means)
  colnames(Hexp_means_df_04) <- c("value", "replic", "samp_size")
  Hexp_means_df_04$samp_size <- 
    factor(Hexp_means_df_04$samp_size, levels = unique(
      as.character(Hexp_means_df_04$samp_size)))
  
  Hexp_means_df_04$marker_num <- "4 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_04[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_04 <- melt(ar_means)
  colnames(ar_means_df_04) <- c("value", "samp_size")
  ar_means_df_04$samp_size <- 
    factor(ar_means_df_04$samp_size, levels = unique(
      as.character(ar_means_df_04$samp_size)))
  
  ar_means_df_04$marker_num <- "4 markers"
  
  # sim_data_03 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_03 <- list()
  for(i in samp_size){ 
    sim_data_03[[i]] <-
      replicate (replic_num, data[[3]][sample(1:nrow(data[[3]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_03[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_03 <- melt(Hobs_means)
  colnames(Hobs_means_df_03) <- c("value", "replic", "samp_size")
  Hobs_means_df_03$samp_size <- 
    factor(Hobs_means_df_03$samp_size, levels = unique(
      as.character(Hobs_means_df_03$samp_size)))
  
  Hobs_means_df_03$marker_num <- "3 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_03 <- melt(Hexp_means)
  colnames(Hexp_means_df_03) <- c("value", "replic", "samp_size")
  Hexp_means_df_03$samp_size <- 
    factor(Hexp_means_df_03$samp_size, levels = unique(
      as.character(Hexp_means_df_03$samp_size)))
  
  Hexp_means_df_03$marker_num <- "3 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_03[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_03 <- melt(ar_means)
  colnames(ar_means_df_03) <- c("value", "samp_size")
  ar_means_df_03$samp_size <- 
    factor(ar_means_df_03$samp_size, levels = unique(
      as.character(ar_means_df_03$samp_size)))
  
  ar_means_df_03$marker_num <- "3 markers"
  
  # sim_data_02 ###############################################
  # Returns a list with as many elements as the different sample size classes.
  # Within those elements reside the replications specified
  sim_data_02 <- list()
  for(i in samp_size){ 
    sim_data_02[[i]] <-
      replicate (replic_num, data[[2]][sample(1:nrow(data[[2]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_02[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot02
  Hobs_means_df_02 <- melt(Hobs_means)
  colnames(Hobs_means_df_02) <- c("value", "replic", "samp_size")
  Hobs_means_df_02$samp_size <- 
    factor(Hobs_means_df_02$samp_size, levels = unique(
      as.character(Hobs_means_df_02$samp_size)))
  
  Hobs_means_df_02$marker_num <- "2 markers"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot02
  Hexp_means_df_02 <- melt(Hexp_means)
  colnames(Hexp_means_df_02) <- c("value", "replic", "samp_size")
  Hexp_means_df_02$samp_size <- 
    factor(Hexp_means_df_02$samp_size, levels = unique(
      as.character(Hexp_means_df_02$samp_size)))

  Hexp_means_df_02$marker_num <- "2 markers"
  
  # Calculation of AR for each generated dataset
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_02[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot02
  ar_means_df_02 <- melt(ar_means)
  colnames(ar_means_df_02) <- c("value", "samp_size")
  ar_means_df_02$samp_size <- 
    factor(ar_means_df_02$samp_size, levels = unique(
      as.character(ar_means_df_02$samp_size)))
  
  ar_means_df_02$marker_num <- "2 markers"
  
  # sim_data_01 ###############################################
  sim_data_01 <- list()
  for(i in samp_size){ 
    sim_data_01[[i]] <-
      replicate (replic_num, data[[1]][sample(1:nrow(data[[1]]$tab), 
                                              i, replace = F)])
  }
  
  # Summaries (Ho, He, allele number) for each generated dataset 
  results <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_data_01[[i]][[j]])
    }
  }
  
  # Ho for each generated dataset 
  Hobs <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hobs[[i]][[j]] <- results[[i]][[j]][["Hobs"]]
    }
  }
  
  # Mean Ho values for each generated dataset 
  Hobs_means <- list()
  for(i in samp_size){
    Hobs_means[[i]] <- lapply(Hobs[[i]], mean)
  }
  
  # Produce a data frame to be plotted by ggplot2
  Hobs_means_df_01 <- melt(Hobs_means)
  colnames(Hobs_means_df_01) <- c("value", "replic", "samp_size")
  Hobs_means_df_01$samp_size <- 
    factor(Hobs_means_df_01$samp_size, levels = unique(
      as.character(Hobs_means_df_01$samp_size)))
  
  Hobs_means_df_01$marker_num <- "1 marker"
  
  # He for each generated dataset 
  Hexp <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      Hexp[[i]][[j]] <- results[[i]][[j]][["Hexp"]]
    }
  }
  
  # Mean He values for each generated dataset 
  Hexp_means <- list()
  for(i in samp_size){
    Hexp_means[[i]] <- lapply(Hexp[[i]], mean)
  }
  
  # Produce a data frame to plotted by ggplot2
  Hexp_means_df_01 <- melt(Hexp_means)
  colnames(Hexp_means_df_01) <- c("value", "replic", "samp_size")
  Hexp_means_df_01$samp_size <- 
    factor(Hexp_means_df_01$samp_size, levels = unique(
      as.character(Hexp_means_df_01$samp_size)))
  
  Hexp_means_df_01$marker_num <- "1 marker"
  
  # Calculation of AR for each generated dataset
  
  # Seperate dataset by pop & select pop to analyze
  obj_fix_list <- seppop(obj_fix) # separate pops
  obj_fix <- obj_fix_list[[pop]] # select pop to analyze
  
  # If there are missing data in a single-locus dataset, poppr package deletes the individual
  # For this reason it might be needed to reduce the highest sample size that can be sampled.
  data_length <- nrow(obj_fix@tab)
  samp_size[length(samp_size)] <- as.character(data_length)
  
  sim_data_01 <- list()
  for(i in samp_size){ 
    sim_data_01[[i]] <-
      replicate (replic_num, obj_fix[sample(1:nrow(obj_fix$tab), 
                                            i, replace = F)])
  }
  
  # Change the name of highest sample size list element and reset samp_size to original value
  real_samp_size <- samp_size[length(samp_size)] 
  if (high_samp_size != real_samp_size){
    sim_data_01[[high_samp_size]] <- sim_data_01[[real_samp_size]]
    sim_data_01[[real_samp_size]] <- NULL
    samp_size[length(samp_size)] <- high_samp_size  # reset to original value
  }
  
  
  ar <- list()
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar[[i]][[j]] <- allelic.richness(sim_data_01[[i]][[j]])
    }
  }
  
  # Mean AR values for each generated dataset
  ar_means <- list() 
  for(i in samp_size){
    for(j in 1:replic_num){ # number of replications
      ar_means[[i]][[j]] <- colMeans(ar[[i]][[j]][["Ar"]])
    }
  }
  
  # Produce a data frame to plotted by ggplot2
  ar_means_df_01 <- melt(ar_means)
  colnames(ar_means_df_01) <- c("value", "samp_size")
  ar_means_df_01$samp_size <- 
    factor(ar_means_df_01$samp_size, levels = unique(
      as.character(ar_means_df_01$samp_size)))
  
  ar_means_df_01$marker_num <- "1 marker"
  
}) 

# Plots ####  

pdf(paste("Abies", pop, "100_repl.pdf", sep = "_"), 
    width = 24, height = 13.5, compress = FALSE)

Hobs_means_tidy <- bind_rows(mget(ls(pattern = "Hobs_means_df_")))

Hobs_means_tidy$marker_num <- 
  factor(Hobs_means_tidy$marker_num, levels = unique(
    as.character(Hobs_means_tidy$marker_num)))

if(pop == "DE_Adult"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (pop == "GR_Adult"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Adult"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (pop == "DE_Regen"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (pop == "GR_Regen"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Regen"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (pop == "DE_Seed"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (pop == "GR_Seed"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Seed"){
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
}else{
  title_Ho <- expression(paste(
    "Observed Heterozygosity (Ho) by sample size & marker number"))}
title_Ho <- expression(paste(
  "Observed Heterozygosity (Ho) by sample size & marker number for german regeneration population of ", 
  italic("A. alba")))

y_axis_Ho <- seq(0.05, 0.95, 0.05)

p_Ho_tidy <- ggplot(Hobs_means_tidy, aes(x = samp_size, y = value)) + 
  geom_boxplot() + 
  facet_wrap(~ marker_num, nrow = 2)

p_Ho_tidy + ggtitle(title_Ho) + xlab("Sample Size") + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Observed Heterozygosity (Ho)", breaks = y_axis_Ho)



Hexp_means_tidy <- bind_rows(mget(ls(pattern = "Hexp_means_df_")))

Hexp_means_tidy$marker_num <- 
  factor(Hexp_means_tidy$marker_num, levels = unique(
    as.character(Hexp_means_tidy$marker_num)))

if(pop == "DE_Adult"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (pop == "GR_Adult"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Adult"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (pop == "DE_Regen"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (pop == "GR_Regen"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Regen"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (pop == "DE_Seed"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (pop == "GR_Seed"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Seed"){
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
}else{
  title_He <- expression(paste(
    "Expected Heterozygosity (He) by sample size & marker number"))}

y_axis_He <- seq(0.05, 0.95, 0.05)

p_He_tidy <- ggplot(Hexp_means_tidy, aes(x = samp_size, y = value)) + 
  geom_boxplot() + 
  facet_wrap(~ marker_num, nrow = 2)

p_He_tidy + ggtitle(title_He) + xlab("Sample Size") + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Expected Heterozygosity (He)", breaks = y_axis_He)



ar_means_tidy <- bind_rows(mget(ls(pattern = "ar_means_df_")))

ar_means_tidy$marker_num <- 
  factor(ar_means_tidy$marker_num, levels = unique(
    as.character(ar_means_tidy$marker_num)))

if(pop == "DE_Adult"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (pop == "GR_Adult"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Adult"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (pop == "DE_Regen"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (pop == "GR_Regen"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Regen"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (pop == "DE_Seed"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (pop == "GR_Seed"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (pop == "SL_Seed"){
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
}else{
  title_ar <- expression(paste(
    "Allelic richness (Ar) by sample size & marker number"))}

y_axis_ar <- 1:30

p_ar_tidy <- ggplot(ar_means_tidy, aes(x = samp_size, y = value)) + 
  geom_boxplot() + 
  facet_wrap(~ marker_num, nrow = 2)

p_ar_tidy + ggtitle(title_ar) + xlab("Sample Size") + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Allelic richness (Ar)", breaks = y_axis_ar)

dev.off()

# Reproducibility ####
writeLines(capture.output(sessionInfo()), paste("sessionInfo_Abies", pop, ".txt", sep = "_"))
