# A GenAlEx formatted excel sheet (.xlsx) is the required input file.

# Set working directory, where the input file is. On RStudio you can achieve 
# this by navigating to "Session" -> "Set Working Directory" -> "Choose Directory"

# Warning of entirely non-type individuals when importing datasets

# After running some of the commands the following warnings appear: 
# "In validityMethod(object) :
# @tab does not contain integers; as of adegenet_2.0-0, numeric values are no longer used"
# They are triggered because of replacing missing data with the mean value. 
# In this case they are harmless.
# https://groups.google.com/forum/#!topic/poppr/F-HImtnMrA8

# Install (if needed) & load required libraries
if (!require("adegenet")) install.packages("adegenet")
library(adegenet)
if (!require("popprxl")) install.packages("popprxl")
library(popprxl)
if (!require("hierfstat")) install.packages("hierfstat")
library(hierfstat)
if (!require("mmod")) install.packages("mmod")
library(mmod)
if (!require("reshape2")) install.packages("reshape2")
library(reshape2)
if (!require("ggplot2")) install.packages("ggplot2") 
library(ggplot2)
if (!require("dplyr")) install.packages("dplyr") 
library(dplyr)
if (!require("tidyr")) install.packages("tidyr") 
library(tidyr)
if (!require("RColorBrewer")) install.packages("RColorBrewer") 
library(RColorBrewer)

set.seed(1994)

# Load dataset
# A GenAlEx formatted excel sheet is the required input
obj <- read.genalexcel(
  "LGM_DE_SI_GR_final.xlsx",   # name of excel file
  sheet = "Abies",             # name of sheet where the genotypes reside
  genclone = F) 

# Real-world datasets are expected to contain missing data which might introduce bias in 
# simulations. With the following function, missing data are replaced with mean values.
# Comment next line to leave missing data as they are.
obj <- missingno(obj, type = "mean")

species <- "Abies"  # "Abies" or "Fagus"
pop <- "GR_Adult" # select pop to analyze (acceptable names: "SL_Adult", "SL_Regen",)

replic_num <- 100   # set number of replications

# Simulations ####
system.time({
  
  # If you'de like to calculate only for EST-SSRs or nSSRs of the Abies dataset uncomment the 
  # appropriate of the following three lines.
  #nSSRs <- c("SF1", "NFF3", "Aag01", "NFH15", "NFF7", "SFb4")
  #EST_SSRs <- c("Aat06", "Aat11", "Aat15", "Aat01", "Aat04")
  #obj <- obj[loc = EST_SSRs]   # EST_SSRs or nSSRs
  
  loci <- sort(nAll(obj)) # vector containing loci from least to
  # most polymorphic according to the pooled dataset from all countries and pops
  loci <- names(loci) # transform named vector to vector of names
  most_poly_locus <- loci[length(loci)]
  
  # Seperate dataset by pop & select pop to analyze
  obj_list <- seppop(obj) # separate pops
  obj <- obj_list[[pop]] 
  id <- paste(species, pop, sep = "_")
  
  
  # There is a bug in the current versions of hierfstat (on cran & development version)
  # which prevents the package from calculating Ar for single-locus datasets.
  # A bug report has been sent: https://github.com/jgx65/hierfstat/issues/25
  # The workaround implemented here, is to create a new dataset where the single locus
  # is duplicated, so mean measures of Ar are unaffected.
  
  # Convert obj to data.frame in order to manipulate it easily
  obj_fix <- genind2df(obj)
  obj_fix <- obj_fix[,c("pop", most_poly_locus)]
  
  # Duplicate column - Create genind object
  duplicate_col <- obj_fix[,2]
  obj_fix$duplicate <- duplicate_col
  
  # Remove pop columns as they are wrongly made into alleles in df2genind
  obj_fix <- obj_fix[-1]
  
  # Create genind object
  obj_fix <- df2genind(obj_fix, sep=NULL, ncode = 3)
  pop(obj_fix) <- rep(pop, nrow(obj_fix@tab))
  
  
  # Set sample size
  if(id %in% c("Abies_DE_Adult", "Abies_GR_Adult", "Fagus_DE_Adult", "Fagus_GR_Adult")){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250")
  }else if (id == "Abies_SL_Adult"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","249")
  }else if (id == "Fagus_SL_Adult"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","251")
  }else if (id %in% c("Abies_DE_Regen", "Abies_GR_Regen", "Abies_SL_Regen",
                      "Fagus_DE_Regen", "Fagus_GR_Regen", "Fagus_SL_Regen")){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200")
  }else if (id %in% c("Abies_DE_Seed","Abies_SL_Seed",
                      "Fagus_DE_Seed", "Fagus_GR_Seed", "Fagus_SL_Seed")){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250","275", "300", "325", "350", "375", "400")
  }else if (id == "Abies_GR_Seed"){
    samp_size <- c("10", "15", "25", "30", "50", "75", "100", "125", "150", "175", "200", 
                   "225","250","275", "300", "325", "350", "375", "382")
  }else{
    print("Unknown population - Cannot continue")
  }
  
  data <- list()
  for(i in 1:length(loci)){
    data[[length(loci)+1-i]] <- obj[, loc = loci[i:length(loci)]]
  }
  
  # Save the highest sample size 
  high_samp_size <- samp_size[length(samp_size)]# added because of hierfstat bug 
  
  # Calculates allelic frequencies for empirical dataset
  f1 <- 0.01 # insert frequency
  f2 <- 0.05
  
  # functions
  sim_dataset_fun <- function(input){
    sim_data <- list()
    for(i in samp_size){ 
      sim_data[[i]] <-
        replicate (replic_num, input[sample(1:nrow(input$tab), 
                                            i, replace = F)])
    }
    return(sim_data)
  }
  
  
  results_fun <- function(sim_dataset){
    results <- list()
    for(i in samp_size){
      for(j in 1:replic_num){ # number of replications
        results[[i]][[j]] <- summary(sim_dataset[[i]][[j]])
      }
    }
    return(results)
  }
  
  Hobs_fun <- function(sim_dataset, results){
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
    Hobs_means_df <- melt(Hobs_means)
    colnames(Hobs_means_df) <- c("value", "replic", "samp_size")
    Hobs_means_df$samp_size <- 
      factor(Hobs_means_df$samp_size, levels = unique(
        as.character(Hobs_means_df$samp_size)))
    
    Hobs_means_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "marker", sep = " ")
    }
    
    return(Hobs_means_df)
  }
  
  
  Hexp_fun <- function(sim_dataset, results){
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
    Hexp_means_df <- melt(Hexp_means)
    colnames(Hexp_means_df) <- c("value", "replic", "samp_size")
    Hexp_means_df$samp_size <- 
      factor(Hexp_means_df$samp_size, levels = unique(
        as.character(Hexp_means_df$samp_size)))
    
    Hexp_means_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "marker", sep = " ")
    }
    
    
    return(Hexp_means_df)
  }
  
  
  ar_fun <- function(sim_dataset){
    # Calculation of AR for each generated dataset
    ar <- list()
    for(i in samp_size){
      for(j in 1:replic_num){ # number of replications
        ar[[i]][[j]] <- allelic.richness(sim_dataset[[i]][[j]])
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
    ar_means_df <- melt(ar_means)
    colnames(ar_means_df) <- c("value", "samp_size")
    ar_means_df$samp_size <- 
      factor(ar_means_df$samp_size, levels = unique(
        as.character(ar_means_df$samp_size)))
    
    ar_means_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "marker", sep = " ")
    }
    
    return(ar_means_df)
  }
  
  
  ar_single_locus_fun <- function(){
    # Calculation of AR for each generated dataset
    
    # Seperate dataset by pop & select pop to analyze
    # obj_fix_list <- seppop(obj_fix) # separate pops
    # obj_fix <- obj_fix_list[[pop]] # select pop to analyze
    
    # If there are missing data in a single-locus dataset, poppr package deletes the individual
    # For this reason it might be needed to reduce the highest sample size that can be sampled.
    data_length <- nrow(obj_fix@tab)
    samp_size[length(samp_size)] <- as.character(data_length)
    
    sim_data_01_fix <- list()
    for(i in samp_size){ 
      sim_data_01_fix[[i]] <-
        replicate (replic_num, obj_fix[sample(1:nrow(obj_fix$tab), 
                                              i, replace = F)])
    }
    
    # Change the name of highest sample size list element and reset samp_size to original value
    real_samp_size <- samp_size[length(samp_size)] 
    if (high_samp_size != real_samp_size){
      sim_data_01_fix[[high_samp_size]] <- sim_data_01_fix[[real_samp_size]]
      sim_data_01_fix[[real_samp_size]] <- NULL
      samp_size[length(samp_size)] <- high_samp_size  # reset to original value
    }
    
    
    ar <- list()
    for(i in samp_size){
      for(j in 1:replic_num){ # number of replications
        ar[[i]][[j]] <- allelic.richness(sim_data_01_fix[[i]][[j]])
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
    
    return(ar_means_df_01)
  }
  
  
  perc_detect_fun <- function(sim_dataset, f1, f2){
    freq_empirical <- apply(tab(sim_dataset[[length(sim_dataset)]][[replic_num]], 
                                freq=TRUE),2,mean, na.rm=TRUE)
    freq_empirical_over_f1 <- subset(freq_empirical, freq_empirical > f1)
    alleles_over_f1 <- names(freq_empirical_over_f1)
    
    freq_empirical_over_f2 <- subset(freq_empirical, freq_empirical > f2)
    alleles_over_f2 <- names(freq_empirical_over_f2)
    
    freq_repl <- list()
    for(i in samp_size){
      for(j in 1:replic_num){
        freq_repl[[i]][[j]] <- apply(tab(sim_dataset[[i]][[j]], 
                                         freq=TRUE),2,mean, na.rm=TRUE)
      }
    }
    
    # Percent of datasets with all alleles f > f1 are found 
    freq_repl_over_f1 <- list()
    freq_over_f1_applied_to_repl <- list()
    n_alleles_over_f1 <- list()
    for(i in samp_size){
      for(j in 1:replic_num){
        freq_repl_over_f1[[i]][[j]] <- freq_repl[[i]][[j]][alleles_over_f1]
        
        freq_over_f1_applied_to_repl[[i]][[j]] <- subset(freq_repl_over_f1[[i]][[j]], 
                                                         freq_repl_over_f1[[i]][[j]] > 0)
        
        n_alleles_over_f1[[i]][[j]] <- length(freq_over_f1_applied_to_repl[[i]][[j]])
      }
    }
    
    for(i in samp_size){
      for(j in 1:replic_num){
        if(n_alleles_over_f1[[i]][[j]] != length(alleles_over_f1)){
          n_alleles_over_f1[[i]][[j]] <- 0
        }else{
          n_alleles_over_f1[[i]][[j]] <- 1
        }
      }
    }
    
    repl_over_f1_detect <- unlist(lapply(n_alleles_over_f1, sum))
    perc_repl_over_f1_detect <- (repl_over_f1_detect / replic_num * 100)
    perc_repl_over_f1_detect <- data.frame("samp_size" = names(perc_repl_over_f1_detect),
                                           "percent_f1" = perc_repl_over_f1_detect)
    perc_repl_over_f1_detect$samp_size <- 
      factor(perc_repl_over_f1_detect$samp_size, levels = unique(
        as.character(perc_repl_over_f1_detect$samp_size)))
    
    
    # Percent of datasets with all alleles f > f2 are found 
    freq_repl_over_f2 <- list()
    freq_over_f2_applied_to_repl <- list()
    n_alleles_over_f2 <- list()
    for(i in samp_size){
      for(j in 1:replic_num){
        freq_repl_over_f2[[i]][[j]] <- freq_repl[[i]][[j]][alleles_over_f2]
        
        freq_over_f2_applied_to_repl[[i]][[j]] <- subset(freq_repl_over_f2[[i]][[j]], 
                                                         freq_repl_over_f2[[i]][[j]] > 0)
        
        n_alleles_over_f2[[i]][[j]] <- length(freq_over_f2_applied_to_repl[[i]][[j]])
      }
    }
    
    for(i in samp_size){
      for(j in 1:replic_num){
        if(n_alleles_over_f2[[i]][[j]] != length(alleles_over_f2)){
          n_alleles_over_f2[[i]][[j]] <- 0
        }else{
          n_alleles_over_f2[[i]][[j]] <- 1
        }
      }
    }
    
    repl_over_f2_detect <- unlist(lapply(n_alleles_over_f2, sum))
    perc_repl_over_f2_detect <- (repl_over_f2_detect / replic_num * 100)
    perc_repl_over_f2_detect <- data.frame("samp_size" = names(perc_repl_over_f2_detect),
                                           "percent_f2" = perc_repl_over_f2_detect)
    perc_repl_over_f2_detect$samp_size <- 
      factor(perc_repl_over_f2_detect$samp_size, levels = unique(
        as.character(perc_repl_over_f2_detect$samp_size)))
    
    perc_repl_detect <- left_join(perc_repl_over_f1_detect, perc_repl_over_f2_detect)
    perc_repl_detect$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "marker", sep = " ")
    }
    
    return(perc_repl_detect)
  }
  
  
  genet_dist_pairs <- function(sim_dataset, empirical, method){
    
    # Rename pop of empirical dataset, so that it can be distinguished from the replicate.
    pop(empirical) <- rep("emp", nrow(empirical@tab))
    
    # Create list with one replicate & the empirical dataset
    pop_pairs_list <- list()
    for(i in samp_size){
      for(j in 1:replic_num){
        pop_pairs_list[[i]][[j]] <- repool(sim_dataset[[i]][[j]], empirical, list = TRUE)
      }
    }
    
    # Create genind objects with one replicate & the empirical dataset
    pop_pairs <- list()
    for(i in samp_size){
      pop_pairs[[i]] <- lapply(pop_pairs_list[[i]], repool)
    }
    
    # Calculate distance
    distance <- list()
    for(i in samp_size){
      distance[[i]] <- as.data.frame(sapply(pop_pairs[[i]], genet.dist, method = method))
    }
    
    for(i in samp_size){
      distance[[i]]$samp_size <- paste(i)
      colnames(distance[[i]]) <- c(method, "samp_size")
    }
    
    # Create a single data.frame with all the distances
    distance_df <- as.data.frame(bind_rows(distance))
    distance_df$samp_size <- factor(distance_df$samp_size, levels = unique(
      as.character(distance_df$samp_size)))
    distance_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])), 
            "marker", sep = " ")
    }
    
    return(distance_df)
  }
  
  
  jost_D_pairs <- function(sim_dataset, empirical){
    
    # Rename pop of empirical dataset, so that it can be distinguished from the replicate.
    pop(empirical) <- rep("emp", nrow(empirical@tab))
    
    # Create list with one replicate & the empirical dataset
    pop_pairs_list <- list()
    for(i in samp_size){
      for(j in 1:replic_num){
        pop_pairs_list[[i]][[j]] <- repool(sim_dataset[[i]][[j]], empirical, list = TRUE)
      }
    }
    
    # Create genind objects with one replicate & the empirical dataset
    pop_pairs <- list()
    for(i in samp_size){
      pop_pairs[[i]] <- lapply(pop_pairs_list[[i]], repool)
    }
    
    # Calculate distance
    distance <- list()
    for(i in samp_size){
      distance[[i]] <- lapply(pop_pairs[[i]], D_Jost)
    }
    
    for(i in samp_size){
      for(j in 1:replic_num){
        distance[[i]][[j]] <- as.data.frame(distance[[i]][[j]][["global.het"]])
      }
    }
    
    for(i in samp_size){
      distance[[i]] <- as.data.frame(bind_rows(distance[[i]]))
      distance[[i]]$samp_size <- paste(i)
    }
    
    
    # Create a single data.frame with all the distances
    distance_df <- as.data.frame(bind_rows(distance))
    colnames(distance_df) <- c("D_Jost", "samp_size")
    distance_df$samp_size <- factor(distance_df$samp_size, levels = unique(
      as.character(distance_df$samp_size)))
    distance_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){
      paste(length(nAll(sim_dataset[[1]][[1]])),
            "markers", sep = " ")
    }else {
      paste(length(nAll(sim_dataset[[1]][[1]])),
            "marker", sep = " ")
    }
    
    return(distance_df)
  }
  
  
  # Simulations
  if(id %in% c("Abies_DE_Adult", "Abies_DE_Regen", "Abies_DE_Seed",
               "Abies_GR_Adult", "Abies_GR_Regen", "Abies_GR_Seed",
               "Fagus_DE_Adult", "Fagus_DE_Regen", "Fagus_DE_Seed")){
    
    sim_data_11 <- sim_dataset_fun(data[[11]])
    # results_11 <- results_fun(sim_data_11)
    # Hobs_means_df_11 <- Hobs_fun(sim_data_11, results_11)
    # Hexp_means_df_11 <- Hexp_fun(sim_data_11, results_11)
    # ar_means_df_11 <- ar_fun(sim_data_11)
    # perc_repl_detect_11 <- perc_detect_fun(sim_data_11, f1, f2)
    fst_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "WC84")
    dch_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "Dch")
    jost_11 <- jost_D_pairs(sim_data_11, data[[11]])
    
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    # results_10 <- results_fun(sim_data_10)
    # Hobs_means_df_10 <- Hobs_fun(sim_data_10, results_10)
    # Hexp_means_df_10 <- Hexp_fun(sim_data_10, results_10)
    # ar_means_df_10 <- ar_fun(sim_data_10)
    # perc_repl_detect_10 <- perc_detect_fun(sim_data_10, f1, f2)
    fst_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "WC84")
    dch_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "Dch")
    jost_10 <- jost_D_pairs(sim_data_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[9]])
    # results_09 <- results_fun(sim_data_09)
    # Hobs_means_df_09 <- Hobs_fun(sim_data_09, results_09)
    # Hexp_means_df_09 <- Hexp_fun(sim_data_09, results_09)
    # ar_means_df_09 <- ar_fun(sim_data_09)
    # perc_repl_detect_09 <- perc_detect_fun(sim_data_09, f1, f2)
    fst_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "WC84")
    dch_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "Dch")
    jost_09 <- jost_D_pairs(sim_data_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[8]])
    # results_08 <- results_fun(sim_data_08)
    # Hobs_means_df_08 <- Hobs_fun(sim_data_08, results_08)
    # Hexp_means_df_08 <- Hexp_fun(sim_data_08, results_08)
    # ar_means_df_08 <- ar_fun(sim_data_08)
    # perc_repl_detect_08 <- perc_detect_fun(sim_data_08, f1, f2)
    fst_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "WC84")
    dch_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "Dch")
    jost_08 <- jost_D_pairs(sim_data_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[7]])
    # results_07 <- results_fun(sim_data_07)
    # Hobs_means_df_07 <- Hobs_fun(sim_data_07, results_07)
    # Hexp_means_df_07 <- Hexp_fun(sim_data_07, results_07)
    # ar_means_df_07 <- ar_fun(sim_data_07)
    # perc_repl_detect_07 <- perc_detect_fun(sim_data_07, f1, f2)
    fst_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "WC84")
    dch_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "Dch")
    jost_07 <- jost_D_pairs(sim_data_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[6]])
    # results_06 <- results_fun(sim_data_06)
    # Hobs_means_df_06 <- Hobs_fun(sim_data_06, results_06)
    # Hexp_means_df_06 <- Hexp_fun(sim_data_06, results_06)
    # ar_means_df_06 <- ar_fun(sim_data_06)
    # perc_repl_detect_06 <- perc_detect_fun(sim_data_06, f1, f2)
    fst_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "WC84")
    dch_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "Dch")
    jost_06 <- jost_D_pairs(sim_data_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[5]])
    # results_05 <- results_fun(sim_data_05)
    # Hobs_means_df_05 <- Hobs_fun(sim_data_05, results_05)
    # Hexp_means_df_05 <- Hexp_fun(sim_data_05, results_05)
    # ar_means_df_05 <- ar_fun(sim_data_05)
    # perc_repl_detect_05 <- perc_detect_fun(sim_data_05, f1, f2)
    fst_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "WC84")
    dch_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "Dch")
    jost_05 <- jost_D_pairs(sim_data_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[4]])
    # results_04 <- results_fun(sim_data_04)
    # Hobs_means_df_04 <- Hobs_fun(sim_data_04, results_04)
    # Hexp_means_df_04 <- Hexp_fun(sim_data_04, results_04)
    # ar_means_df_04 <- ar_fun(sim_data_04)
    # perc_repl_detect_04 <- perc_detect_fun(sim_data_04, f1, f2)
    fst_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "WC84")
    dch_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "Dch")
    jost_04 <- jost_D_pairs(sim_data_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[3]])
    # results_03 <- results_fun(sim_data_03)
    # Hobs_means_df_03 <- Hobs_fun(sim_data_03, results_03)
    # Hexp_means_df_03 <- Hexp_fun(sim_data_03, results_03)
    # ar_means_df_03 <- ar_fun(sim_data_03)
    # perc_repl_detect_03 <- perc_detect_fun(sim_data_03, f1, f2)
    fst_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "WC84")
    dch_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "Dch")
    jost_03 <- jost_D_pairs(sim_data_03, data[[03]])
    
    sim_data_02 <- sim_dataset_fun(data[[2]])
    # results_02 <- results_fun(sim_data_02)
    # Hobs_means_df_02 <- Hobs_fun(sim_data_02, results_02)
    # Hexp_means_df_02 <- Hexp_fun(sim_data_02, results_02)
    # ar_means_df_02 <- ar_fun(sim_data_02)
    # perc_repl_detect_02 <- perc_detect_fun(sim_data_02, f1, f2)
    fst_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "WC84")
    dch_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "Dch")
    jost_02 <- jost_D_pairs(sim_data_02, data[[02]])
    
    # sim_data_01 <- sim_dataset_fun(data[[1]])
    # results_01 <- results_fun(sim_data_01)
    # Hobs_means_df_01 <- Hobs_fun(sim_data_01, results_01)
    # Hexp_means_df_01 <- Hexp_fun(sim_data_01, results_01)
    # perc_repl_detect_01 <- perc_detect_fun(sim_data_01, f1, f2)
    # ar_means_df_01 <- ar_single_locus_fun()
    # fst_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "WC84")
    # dch_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "Dch")
    # jost_01 <- d_jost_pairs(sim_data_01, data[[01]])
    
    
  }else if(id %in% c("Abies_SL_Adult", "Abies_SL_Regen", "Abies_SL_Seed")){
    
    sim_data_17 <- sim_dataset_fun(data[[17]])
    # results_17 <- results_fun(sim_data_17)
    # Hobs_means_df_17 <- Hobs_fun(sim_data_17, results_17)
    # Hexp_means_df_17 <- Hexp_fun(sim_data_17, results_17)
    # ar_means_df_17 <- ar_fun(sim_data_17)
    # perc_repl_detect_17 <- perc_detect_fun(sim_data_17, f1, f2)
    fst_17 <- genet_dist_pairs(sim_data_17, data[[17]], method = "WC84")
    dch_17 <- genet_dist_pairs(sim_data_17, data[[17]], method = "Dch")
    jost_17 <- jost_D_pairs(sim_data_17, data[[17]])
    
    sim_data_16 <- sim_dataset_fun(data[[16]])
    # results_16 <- results_fun(sim_data_16)
    # Hobs_means_df_16 <- Hobs_fun(sim_data_16, results_16)
    # Hexp_means_df_16 <- Hexp_fun(sim_data_16, results_16)
    # ar_means_df_16 <- ar_fun(sim_data_16)
    # perc_repl_detect_16 <- perc_detect_fun(sim_data_16, f1, f2)
    fst_16 <- genet_dist_pairs(sim_data_16, data[[16]], method = "WC84")
    dch_16 <- genet_dist_pairs(sim_data_16, data[[16]], method = "Dch")
    jost_16 <- jost_D_pairs(sim_data_16, data[[16]])
    
    sim_data_15 <- sim_dataset_fun(data[[15]])
    # results_15 <- results_fun(sim_data_15)
    # Hobs_means_df_15 <- Hobs_fun(sim_data_15, results_15)
    # Hexp_means_df_15 <- Hexp_fun(sim_data_15, results_15)
    # ar_means_df_15 <- ar_fun(sim_data_15)
    # perc_repl_detect_15 <- perc_detect_fun(sim_data_15, f1, f2)
    fst_15 <- genet_dist_pairs(sim_data_15, data[[15]], method = "WC84")
    dch_15 <- genet_dist_pairs(sim_data_15, data[[15]], method = "Dch")
    jost_15 <- jost_D_pairs(sim_data_15, data[[15]])
    
    sim_data_14 <- sim_dataset_fun(data[[14]])
    # results_14 <- results_fun(sim_data_14)
    # Hobs_means_df_14 <- Hobs_fun(sim_data_14, results_14)
    # Hexp_means_df_14 <- Hexp_fun(sim_data_14, results_14)
    # ar_means_df_14 <- ar_fun(sim_data_14)
    # perc_repl_detect_14 <- perc_detect_fun(sim_data_14, f1, f2)
    fst_14 <- genet_dist_pairs(sim_data_14, data[[14]], method = "WC84")
    dch_14 <- genet_dist_pairs(sim_data_14, data[[14]], method = "Dch")
    jost_14 <- jost_D_pairs(sim_data_14, data[[14]])
    
    sim_data_13 <- sim_dataset_fun(data[[13]])
    # results_13 <- results_fun(sim_data_13)
    # Hobs_means_df_13 <- Hobs_fun(sim_data_13, results_13)
    # Hexp_means_df_13 <- Hexp_fun(sim_data_13, results_13)
    # ar_means_df_13 <- ar_fun(sim_data_13)
    # perc_repl_detect_13 <- perc_detect_fun(sim_data_13, f1, f2)
    fst_13 <- genet_dist_pairs(sim_data_13, data[[13]], method = "WC84")
    dch_13 <- genet_dist_pairs(sim_data_13, data[[13]], method = "Dch")
    jost_13 <- jost_D_pairs(sim_data_13, data[[13]])
    
    sim_data_12 <- sim_dataset_fun(data[[12]])
    # results_12 <- results_fun(sim_data_12)
    # Hobs_means_df_12 <- Hobs_fun(sim_data_12, results_12)
    # Hexp_means_df_12 <- Hexp_fun(sim_data_12, results_12)
    # ar_means_df_12 <- ar_fun(sim_data_12)
    # perc_repl_detect_12 <- perc_detect_fun(sim_data_12, f1, f2)
    fst_12 <- genet_dist_pairs(sim_data_12, data[[12]], method = "WC84")
    dch_12 <- genet_dist_pairs(sim_data_12, data[[12]], method = "Dch")
    jost_12 <- jost_D_pairs(sim_data_12, data[[12]])
    
    sim_data_11 <- sim_dataset_fun(data[[11]])
    # results_11 <- results_fun(sim_data_11)
    # Hobs_means_df_11 <- Hobs_fun(sim_data_11, results_11)
    # Hexp_means_df_11 <- Hexp_fun(sim_data_11, results_11)
    # ar_means_df_11 <- ar_fun(sim_data_11)
    # perc_repl_detect_11 <- perc_detect_fun(sim_data_11, f1, f2)
    fst_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "WC84")
    dch_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "Dch")
    jost_11 <- jost_D_pairs(sim_data_11, data[[11]])
    
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    # results_10 <- results_fun(sim_data_10)
    # Hobs_means_df_10 <- Hobs_fun(sim_data_10, results_10)
    # Hexp_means_df_10 <- Hexp_fun(sim_data_10, results_10)
    # ar_means_df_10 <- ar_fun(sim_data_10)
    # perc_repl_detect_10 <- perc_detect_fun(sim_data_10, f1, f2)
    fst_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "WC84")
    dch_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "Dch")
    jost_10 <- jost_D_pairs(sim_data_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[9]])
    # results_09 <- results_fun(sim_data_09)
    # Hobs_means_df_09 <- Hobs_fun(sim_data_09, results_09)
    # Hexp_means_df_09 <- Hexp_fun(sim_data_09, results_09)
    # ar_means_df_09 <- ar_fun(sim_data_09)
    # perc_repl_detect_09 <- perc_detect_fun(sim_data_09, f1, f2)
    fst_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "WC84")
    dch_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "Dch")
    jost_09 <- jost_D_pairs(sim_data_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[8]])
    # results_08 <- results_fun(sim_data_08)
    # Hobs_means_df_08 <- Hobs_fun(sim_data_08, results_08)
    # Hexp_means_df_08 <- Hexp_fun(sim_data_08, results_08)
    # ar_means_df_08 <- ar_fun(sim_data_08)
    # perc_repl_detect_08 <- perc_detect_fun(sim_data_08, f1, f2)
    fst_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "WC84")
    dch_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "Dch")
    jost_08 <- jost_D_pairs(sim_data_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[7]])
    # results_07 <- results_fun(sim_data_07)
    # Hobs_means_df_07 <- Hobs_fun(sim_data_07, results_07)
    # Hexp_means_df_07 <- Hexp_fun(sim_data_07, results_07)
    # ar_means_df_07 <- ar_fun(sim_data_07)
    # perc_repl_detect_07 <- perc_detect_fun(sim_data_07, f1, f2)
    fst_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "WC84")
    dch_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "Dch")
    jost_07 <- jost_D_pairs(sim_data_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[6]])
    # results_06 <- results_fun(sim_data_06)
    # Hobs_means_df_06 <- Hobs_fun(sim_data_06, results_06)
    # Hexp_means_df_06 <- Hexp_fun(sim_data_06, results_06)
    # ar_means_df_06 <- ar_fun(sim_data_06)
    # perc_repl_detect_06 <- perc_detect_fun(sim_data_06, f1, f2)
    fst_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "WC84")
    dch_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "Dch")
    jost_06 <- jost_D_pairs(sim_data_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[5]])
    # results_05 <- results_fun(sim_data_05)
    # Hobs_means_df_05 <- Hobs_fun(sim_data_05, results_05)
    # Hexp_means_df_05 <- Hexp_fun(sim_data_05, results_05)
    # ar_means_df_05 <- ar_fun(sim_data_05)
    # perc_repl_detect_05 <- perc_detect_fun(sim_data_05, f1, f2)
    fst_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "WC84")
    dch_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "Dch")
    jost_05 <- jost_D_pairs(sim_data_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[4]])
    # results_04 <- results_fun(sim_data_04)
    # Hobs_means_df_04 <- Hobs_fun(sim_data_04, results_04)
    # Hexp_means_df_04 <- Hexp_fun(sim_data_04, results_04)
    # ar_means_df_04 <- ar_fun(sim_data_04)
    # perc_repl_detect_04 <- perc_detect_fun(sim_data_04, f1, f2)
    fst_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "WC84")
    dch_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "Dch")
    jost_04 <- jost_D_pairs(sim_data_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[3]])
    # results_03 <- results_fun(sim_data_03)
    # Hobs_means_df_03 <- Hobs_fun(sim_data_03, results_03)
    # Hexp_means_df_03 <- Hexp_fun(sim_data_03, results_03)
    # ar_means_df_03 <- ar_fun(sim_data_03)
    # perc_repl_detect_03 <- perc_detect_fun(sim_data_03, f1, f2)
    fst_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "WC84")
    dch_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "Dch")
    jost_03 <- jost_D_pairs(sim_data_03, data[[03]])
    
    sim_data_02 <- sim_dataset_fun(data[[2]])
    # results_02 <- results_fun(sim_data_02)
    # Hobs_means_df_02 <- Hobs_fun(sim_data_02, results_02)
    # Hexp_means_df_02 <- Hexp_fun(sim_data_02, results_02)
    # ar_means_df_02 <- ar_fun(sim_data_02)
    # perc_repl_detect_02 <- perc_detect_fun(sim_data_02, f1, f2)
    fst_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "WC84")
    dch_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "Dch")
    jost_02 <- jost_D_pairs(sim_data_02, data[[02]])
    
    # sim_data_01 <- sim_dataset_fun(data[[1]])
    # results_01 <- results_fun(sim_data_01)
    # Hobs_means_df_01 <- Hobs_fun(sim_data_01, results_01)
    # Hexp_means_df_01 <- Hexp_fun(sim_data_01, results_01)
    # perc_repl_detect_01 <- perc_detect_fun(sim_data_01, f1, f2)
    # ar_means_df_01 <- ar_single_locus_fun()
    # fst_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "WC84")
    # dch_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "Dch")
    # jost_01 <- d_jost_pairs(sim_data_01, data[[01]])
    
    
    
    
  }else if(id %in% c("Fagus_GR_Adult", "Fagus_GR_Regen", "Fagus_GR_Seed", 
                     "Fagus_SL_Adult", "Fagus_SL_Regen", "Fagus_SL_Seed")){
    sim_data_11 <- sim_dataset_fun(data[[11]])
    # results_11 <- results_fun(sim_data_11)
    # Hobs_means_df_11 <- Hobs_fun(sim_data_11, results_11)
    # Hexp_means_df_11 <- Hexp_fun(sim_data_11, results_11)
    # ar_means_df_11 <- ar_fun(sim_data_11)
    # perc_repl_detect_11 <- perc_detect_fun(sim_data_11, f1, f2)
    fst_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "WC84")
    dch_11 <- genet_dist_pairs(sim_data_11, data[[11]], method = "Dch")
    jost_11 <- jost_D_pairs(sim_data_11, data[[11]])
    
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    # results_10 <- results_fun(sim_data_10)
    # Hobs_means_df_10 <- Hobs_fun(sim_data_10, results_10)
    # Hexp_means_df_10 <- Hexp_fun(sim_data_10, results_10)
    # ar_means_df_10 <- ar_fun(sim_data_10)
    # perc_repl_detect_10 <- perc_detect_fun(sim_data_10, f1, f2)
    fst_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "WC84")
    dch_10 <- genet_dist_pairs(sim_data_10, data[[10]], method = "Dch")
    jost_10 <- jost_D_pairs(sim_data_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[9]])
    # results_09 <- results_fun(sim_data_09)
    # Hobs_means_df_09 <- Hobs_fun(sim_data_09, results_09)
    # Hexp_means_df_09 <- Hexp_fun(sim_data_09, results_09)
    # ar_means_df_09 <- ar_fun(sim_data_09)
    # perc_repl_detect_09 <- perc_detect_fun(sim_data_09, f1, f2)
    fst_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "WC84")
    dch_09 <- genet_dist_pairs(sim_data_09, data[[09]], method = "Dch")
    jost_09 <- jost_D_pairs(sim_data_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[8]])
    # results_08 <- results_fun(sim_data_08)
    # Hobs_means_df_08 <- Hobs_fun(sim_data_08, results_08)
    # Hexp_means_df_08 <- Hexp_fun(sim_data_08, results_08)
    # ar_means_df_08 <- ar_fun(sim_data_08)
    # perc_repl_detect_08 <- perc_detect_fun(sim_data_08, f1, f2)
    fst_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "WC84")
    dch_08 <- genet_dist_pairs(sim_data_08, data[[08]], method = "Dch")
    jost_08 <- jost_D_pairs(sim_data_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[7]])
    # results_07 <- results_fun(sim_data_07)
    # Hobs_means_df_07 <- Hobs_fun(sim_data_07, results_07)
    # Hexp_means_df_07 <- Hexp_fun(sim_data_07, results_07)
    # ar_means_df_07 <- ar_fun(sim_data_07)
    # perc_repl_detect_07 <- perc_detect_fun(sim_data_07, f1, f2)
    fst_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "WC84")
    dch_07 <- genet_dist_pairs(sim_data_07, data[[07]], method = "Dch")
    jost_07 <- jost_D_pairs(sim_data_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[6]])
    # results_06 <- results_fun(sim_data_06)
    # Hobs_means_df_06 <- Hobs_fun(sim_data_06, results_06)
    # Hexp_means_df_06 <- Hexp_fun(sim_data_06, results_06)
    # ar_means_df_06 <- ar_fun(sim_data_06)
    # perc_repl_detect_06 <- perc_detect_fun(sim_data_06, f1, f2)
    fst_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "WC84")
    dch_06 <- genet_dist_pairs(sim_data_06, data[[06]], method = "Dch")
    jost_06 <- jost_D_pairs(sim_data_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[5]])
    # results_05 <- results_fun(sim_data_05)
    # Hobs_means_df_05 <- Hobs_fun(sim_data_05, results_05)
    # Hexp_means_df_05 <- Hexp_fun(sim_data_05, results_05)
    # ar_means_df_05 <- ar_fun(sim_data_05)
    # perc_repl_detect_05 <- perc_detect_fun(sim_data_05, f1, f2)
    fst_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "WC84")
    dch_05 <- genet_dist_pairs(sim_data_05, data[[05]], method = "Dch")
    jost_05 <- jost_D_pairs(sim_data_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[4]])
    # results_04 <- results_fun(sim_data_04)
    # Hobs_means_df_04 <- Hobs_fun(sim_data_04, results_04)
    # Hexp_means_df_04 <- Hexp_fun(sim_data_04, results_04)
    # ar_means_df_04 <- ar_fun(sim_data_04)
    # perc_repl_detect_04 <- perc_detect_fun(sim_data_04, f1, f2)
    fst_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "WC84")
    dch_04 <- genet_dist_pairs(sim_data_04, data[[04]], method = "Dch")
    jost_04 <- jost_D_pairs(sim_data_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[3]])
    # results_03 <- results_fun(sim_data_03)
    # Hobs_means_df_03 <- Hobs_fun(sim_data_03, results_03)
    # Hexp_means_df_03 <- Hexp_fun(sim_data_03, results_03)
    # ar_means_df_03 <- ar_fun(sim_data_03)
    # perc_repl_detect_03 <- perc_detect_fun(sim_data_03, f1, f2)
    fst_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "WC84")
    dch_03 <- genet_dist_pairs(sim_data_03, data[[03]], method = "Dch")
    jost_03 <- jost_D_pairs(sim_data_03, data[[03]])
    
    sim_data_02 <- sim_dataset_fun(data[[2]])
    # results_02 <- results_fun(sim_data_02)
    # Hobs_means_df_02 <- Hobs_fun(sim_data_02, results_02)
    # Hexp_means_df_02 <- Hexp_fun(sim_data_02, results_02)
    # ar_means_df_02 <- ar_fun(sim_data_02)
    # perc_repl_detect_02 <- perc_detect_fun(sim_data_02, f1, f2)
    fst_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "WC84")
    dch_02 <- genet_dist_pairs(sim_data_02, data[[02]], method = "Dch")
    jost_02 <- jost_D_pairs(sim_data_02, data[[02]])
    
    # sim_data_01 <- sim_dataset_fun(data[[1]])
    # results_01 <- results_fun(sim_data_01)
    # Hobs_means_df_01 <- Hobs_fun(sim_data_01, results_01)
    # Hexp_means_df_01 <- Hexp_fun(sim_data_01, results_01)
    # perc_repl_detect_01 <- perc_detect_fun(sim_data_01, f1, f2)
    # ar_means_df_01 <- ar_single_locus_fun()
    # fst_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "WC84")
    # dch_01 <- genet_dist_pairs(sim_data_01, data[[01]], method = "Dch")
    # jost_01 <- d_jost_pairs(sim_data_01, data[[01]])
    
  }else{
    print("Unknown population - Cannot continue")
  }
}) 

# Plots ####  

pdf(paste(id, "100_repl.pdf", sep = "_"), 
    width = 24, height = 13.5, compress = FALSE)

my_palette <- brewer.pal(12, "Set3") # create a new palette
my_palette <- colorRampPalette(my_palette)(19) # how many colors this palette will have


fst_tidy <- bind_rows(mget(ls(pattern = "fst_")))
fst_tidy$marker_num <- 
  factor(fst_tidy$marker_num, levels = unique(
    as.character(fst_tidy$marker_num)))

if(id == "Abies_DE_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number"))}

y_axis_fst <- seq(0, 0.99, 0.001)

p_fst_tidy <- ggplot(fst_tidy, aes(x = samp_size, y = WC84)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_fst_tidy + ggtitle(title_Fst) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Mean pairwise Fst", breaks = y_axis_fst) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")



dch_tidy <- bind_rows(mget(ls(pattern = "dch_")))
dch_tidy$marker_num <- 
  factor(dch_tidy$marker_num, levels = unique(
    as.character(dch_tidy$marker_num)))

if(id == "Abies_DE_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number"))}

y_axis_dch <- seq(0, 0.99, 0.01)

p_dch_tidy <- ggplot(dch_tidy, aes(x = samp_size, y = Dch)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_dch_tidy + ggtitle(title_dch) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Cavalli-Sforza and Edwards Chord distance", breaks = y_axis_dch) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")



jost_tidy <- bind_rows(mget(ls(pattern = "jost_")))
jost_tidy$marker_num <- 
  factor(jost_tidy$marker_num, levels = unique(
    as.character(jost_tidy$marker_num)))

if(id == "Abies_DE_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number"))}

y_axis_jost <- seq(0, 0.99, 0.01)

p_jost_tidy <- ggplot(jost_tidy, aes(x = samp_size, y = D_Jost)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_jost_tidy + ggtitle(title_jost) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Jost's D", breaks = y_axis_jost) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")

dev.off()

# Reproducibility ####
writeLines(capture.output(sessionInfo()), paste("sessionInfo", id, ".txt", sep = "_"))