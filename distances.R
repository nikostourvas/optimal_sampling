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
if (!require("devtools")) install.packages("devtools") 
library(devtools)
if (!require("hierfstat")) install_github("jgx65/hierfstat") 
library("hierfstat")

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
pop <- "GR_Seed" # select pop to analyze (acceptable names: "SL_Adult", "SL_Regen",)

replic_num <- 100   # set number of replications

# Simulations ####
system.time({

loci <- sort(nAll(obj)) # vector containing loci from least to
# most polymorphic according to the pooled dataset from all countries and pops
loci <- names(loci) # transform named vector to vector of names
most_poly_locus <- loci[length(loci)]

# Seperate dataset by pop & select pop to analyze
obj_list <- seppop(obj) # separate pops
obj <- obj_list[[pop]] 
id <- paste(species, pop, sep = "_")


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

# sim_dataset_fun without creating multiple iterations of the empirical data set
sim_dataset_fun <- function(empirical){
  sim_data <- list()
  for(i in samp_size[-length(samp_size)]){ 
    sim_data[[i]] <-
      replicate (replic_num, empirical[sample(1:nrow(empirical$tab), 
                                          i, replace = F)])
  }
  
  sim_data[[samp_size[length(samp_size)]]] <- empirical
  return(sim_data)
}

pairs_creator <- function(sim_dataset, empirical){
  
  # Rename pop of empirical dataset, so that it can be distinguished from the replicate.
  pop(empirical) <- rep("emp", nrow(empirical@tab))
  
  # Create list with one replicate & the empirical dataset
  pop_pairs_list <- list()
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){
      pop_pairs_list[[i]][[j]] <- repool(sim_dataset[[i]][[j]], empirical, list = TRUE)
    }
  }
  
  pop_pairs_list[[samp_size[length(samp_size)]]] <- 
    repool(sim_dataset[[samp_size[length(samp_size)]]], empirical, list = FALSE)
  
  
  # Create genind objects with one replicate & the empirical dataset
  pop_pairs <- list()
  for(i in samp_size[-length(samp_size)]){
    pop_pairs[[i]] <- lapply(pop_pairs_list[[i]], repool)
  }
  
  pop_pairs[[samp_size[length(samp_size)]]] <- pop_pairs_list[[samp_size[length(samp_size)]]]
  
  return(pop_pairs)
}



genet_dist_pairs <- function(pop_pairs, method){
  
  # Calculate distance
  distance <- list()
  for(i in samp_size[-length(samp_size)]){
    distance[[i]] <- as.data.frame(sapply(pop_pairs[[i]], genet.dist, method = method))
  }
  
  # for empirical data set 
  distance[[samp_size[length(samp_size)]]] <- 
    as.vector(genet.dist(pop_pairs[[samp_size[length(samp_size)]]]))
  
  # for empirical data set 
  distance[[samp_size[length(samp_size)]]]$samp_size <- paste(samp_size[length(samp_size)])
  distance[[samp_size[length(samp_size)]]] <- as.data.frame(distance[[samp_size[length(samp_size)]]])
  
  for(i in samp_size[-length(samp_size)]){
    distance[[i]]$samp_size <- paste(i)
  }
  
  for(i in samp_size){
    colnames(distance[[i]]) <- c(method, "samp_size")
  }
  
  # Create a single data.frame with all the distances
  distance_df <- as.data.frame(bind_rows(distance))
  distance_df$samp_size <- factor(distance_df$samp_size, levels = unique(
    as.character(distance_df$samp_size)))
  distance_df$marker_num <- if(length(nAll(pop_pairs[[1]][[1]])) > 1){ 
    paste(length(nAll(pop_pairs[[1]][[1]])), 
          "markers", sep = " ")
  }else {
    paste(length(nAll(pop_pairs[[1]][[1]])), 
          "marker", sep = " ")
  }
  
  # Replace negative values with zero
  distance_df$original_values <- distance_df[, method] # keep original values
  distance_df[, method][distance_df[, method] < 0] <- 0
  
  
  return(distance_df)
}




jostD_pairs <- function(pop_pairs, empirical){
  
  # Calculate distance
  distance <- list()
  for(i in samp_size[-length(samp_size)]){
    distance[[i]] <- lapply(pop_pairs[[i]], D_Jost)
  }
  
  # for empirical data set
  distance[[samp_size[length(samp_size)]]] <- D_Jost(pop_pairs[[samp_size[length(samp_size)]]])
  
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){
      distance[[i]][[j]] <- as.data.frame(distance[[i]][[j]][["global.het"]])
    }
  }
  
  # for empirical data set
  distance[[samp_size[length(samp_size)]]] <- as.data.frame(
    distance[[samp_size[length(samp_size)]]][["global.het"]])
  
  for(i in samp_size){
    distance[[i]] <- as.data.frame(bind_rows(distance[[i]]))
    distance[[i]]$samp_size <- paste(i)
  }
  
  # for empirical data set
  colnames(distance[[samp_size[length(samp_size)]]]) <- colnames(distance[[1]])
  
  # Create a single data.frame with all the distances
  distance_df <- as.data.frame(bind_rows(distance))
  colnames(distance_df) <- c("D_Jost", "samp_size")
  distance_df$samp_size <- factor(distance_df$samp_size, levels = unique(
    as.character(distance_df$samp_size)))
  distance_df$marker_num <- if(length(nAll(pop_pairs[[1]][[1]])) > 1){
    paste(length(nAll(pop_pairs[[1]][[1]])),
          "markers", sep = " ")
  }else {
    paste(length(nAll(pop_pairs[[1]][[1]])),
          "marker", sep = " ")
  }
  
  # Replace negative values with zero
  distance_df$original_values <- distance_df[, "D_Jost"] # keep original values
  distance_df[, "D_Jost"][distance_df[, "D_Jost"] < 0] <- 0
  
  return(distance_df)
}
  


Gst_pairs <- function(pop_pairs, empirical){
  
  # Calculate distance
  distance <- list()
  for(i in samp_size[-length(samp_size)]){
    distance[[i]] <- lapply(pop_pairs[[i]], Gst_Nei)
  }
  
  # for empirical data set
  distance[[samp_size[length(samp_size)]]] <- Gst_Nei(pop_pairs[[samp_size[length(samp_size)]]])
  
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){
      distance[[i]][[j]] <- as.data.frame(distance[[i]][[j]][["global"]])
    }
  }
  
  # for empirical data set
  distance[[samp_size[length(samp_size)]]] <- as.data.frame(
    distance[[samp_size[length(samp_size)]]][["global"]])
  
  for(i in samp_size){
    distance[[i]] <- as.data.frame(bind_rows(distance[[i]]))
    distance[[i]]$samp_size <- paste(i)
  }
  
  # for empirical data set
  colnames(distance[[samp_size[length(samp_size)]]]) <- colnames(distance[[1]])
  
  # Create a single data.frame with all the distances
  distance_df <- as.data.frame(bind_rows(distance))
  colnames(distance_df) <- c("Gst_Nei", "samp_size")
  distance_df$samp_size <- factor(distance_df$samp_size, levels = unique(
    as.character(distance_df$samp_size)))
  distance_df$marker_num <- if(length(nAll(pop_pairs[[1]][[1]])) > 1){
    paste(length(nAll(pop_pairs[[1]][[1]])),
          "markers", sep = " ")
  }else {
    paste(length(nAll(pop_pairs[[1]][[1]])),
          "marker", sep = " ")
  }
  
  # Replace negative values with zero
  distance_df$original_values <- distance_df[, "Gst_Nei"] # keep original values
  distance_df[, "Gst_Nei"][distance_df[, "Gst_Nei"] < 0] <- 0
  
  return(distance_df)
}
  
  # Simulations
  if(id %in% c("Abies_DE_Adult", "Abies_DE_Regen", "Abies_DE_Seed",
               "Abies_GR_Adult", "Abies_GR_Regen", "Abies_GR_Seed",
               "Fagus_DE_Adult", "Fagus_DE_Regen", "Fagus_DE_Seed")){
    
    sim_data_11 <- sim_dataset_fun(data[[11]])
    pop_pairs_11 <- pairs_creator(sim_data_11, data[[11]])
    fst_11 <- genet_dist_pairs(pop_pairs_11, method = "WC84")
    dch_11 <- genet_dist_pairs(pop_pairs_11, method = "Dch")
    jost_11 <- jostD_pairs(pop_pairs_11, data[[11]])
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    pop_pairs_10 <-pairs_creator(sim_data_10, data[[10]])
    fst_10 <- genet_dist_pairs(pop_pairs_10, method = "WC84")
    dch_10 <- genet_dist_pairs(pop_pairs_10, method = "Dch")
    jost_10 <- jostD_pairs(pop_pairs_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[09]])
    pop_pairs_09 <- pairs_creator(sim_data_09, data[[09]])
    fst_09 <- genet_dist_pairs(pop_pairs_09, method = "WC84")
    dch_09 <- genet_dist_pairs(pop_pairs_09, method = "Dch")
    jost_09 <- jostD_pairs(pop_pairs_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[08]])
    pop_pairs_08 <- pairs_creator(sim_data_08, data[[08]])
    fst_08 <- genet_dist_pairs(pop_pairs_08, method = "WC84")
    dch_08 <- genet_dist_pairs(pop_pairs_08, method = "Dch")
    jost_08 <- jostD_pairs(pop_pairs_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[07]])
    pop_pairs_07 <- pairs_creator(sim_data_07, data[[07]])
    fst_07 <- genet_dist_pairs(pop_pairs_07, method = "WC84")
    dch_07 <- genet_dist_pairs(pop_pairs_07, method = "Dch")
    jost_07 <- jostD_pairs(pop_pairs_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[06]])
    pop_pairs_06 <- pairs_creator(sim_data_06, data[[06]])
    fst_06 <- genet_dist_pairs(pop_pairs_06, method = "WC84")
    dch_06 <- genet_dist_pairs(pop_pairs_06, method = "Dch")
    jost_06 <- jostD_pairs(pop_pairs_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[05]])
    pop_pairs_05 <- pairs_creator(sim_data_05, data[[05]])
    fst_05 <- genet_dist_pairs(pop_pairs_05, method = "WC84")
    dch_05 <- genet_dist_pairs(pop_pairs_05, method = "Dch")
    jost_05 <- jostD_pairs(pop_pairs_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[04]])
    pop_pairs_04 <- pairs_creator(sim_data_04, data[[04]])
    fst_04 <- genet_dist_pairs(pop_pairs_04, method = "WC84")
    dch_04 <- genet_dist_pairs(pop_pairs_04, method = "Dch")
    jost_04 <- jostD_pairs(pop_pairs_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[03]])
    pop_pairs_03 <- pairs_creator(sim_data_03, data[[03]])
    fst_03 <- genet_dist_pairs(pop_pairs_03, method = "WC84")
    dch_03 <- genet_dist_pairs(pop_pairs_03, method = "Dch")
    jost_03 <- jostD_pairs(pop_pairs_03, data[[03]])

    sim_data_02 <- sim_dataset_fun(data[[02]])
    pop_pairs_02 <- pairs_creator(sim_data_02, data[[02]])
    fst_02 <- genet_dist_pairs(pop_pairs_02, method = "WC84")
    dch_02 <- genet_dist_pairs(pop_pairs_02, method = "Dch")
    jost_02 <- jostD_pairs(pop_pairs_02, data[[02]])

    
  }else if(id %in% c("Abies_SL_Adult", "Abies_SL_Regen", "Abies_SL_Seed")){
    
    sim_data_17 <- sim_dataset_fun(data[[17]])
    pop_pairs_17 <- pairs_creator(sim_data_17, data[[17]])
    fst_17 <- genet_dist_pairs(pop_pairs_17, method = "WC84")
    dch_17 <- genet_dist_pairs(pop_pairs_17, method = "Dch")
    jost_17 <- jostD_pairs(pop_pairs_17, data[[17]])
    
    sim_data_16 <- sim_dataset_fun(data[[16]])
    pop_pairs_16 <- pairs_creator(sim_data_16, data[[16]])
    fst_16 <- genet_dist_pairs(pop_pairs_16, method = "WC84")
    dch_16 <- genet_dist_pairs(pop_pairs_16, method = "Dch")
    jost_16 <- jostD_pairs(pop_pairs_16, data[[16]])
    
    sim_data_15 <- sim_dataset_fun(data[[15]])
    pop_pairs_15 <- pairs_creator(sim_data_15, data[[15]])
    fst_15 <- genet_dist_pairs(pop_pairs_15, method = "WC84")
    dch_15 <- genet_dist_pairs(pop_pairs_15, method = "Dch")
    jost_15 <- jostD_pairs(pop_pairs_15, data[[15]])
    
    sim_data_14 <- sim_dataset_fun(data[[14]])
    pop_pairs_14 <- pairs_creator(sim_data_14, data[[14]])
    fst_14 <- genet_dist_pairs(pop_pairs_14, method = "WC84")
    dch_14 <- genet_dist_pairs(pop_pairs_14, method = "Dch")
    jost_14 <- jostD_pairs(pop_pairs_14, data[[14]])
    
    sim_data_13 <- sim_dataset_fun(data[[13]])
    pop_pairs_13 <- pairs_creator(sim_data_13, data[[13]])
    fst_13 <- genet_dist_pairs(pop_pairs_13, method = "WC84")
    dch_13 <- genet_dist_pairs(pop_pairs_13, method = "Dch")
    jost_13 <- jostD_pairs(pop_pairs_13, data[[13]])
    
    sim_data_12 <- sim_dataset_fun(data[[12]])
    pop_pairs_12 <- pairs_creator(sim_data_12, data[[12]])
    fst_12 <- genet_dist_pairs(pop_pairs_12, method = "WC84")
    dch_12 <- genet_dist_pairs(pop_pairs_12, method = "Dch")
    jost_12 <- jostD_pairs(pop_pairs_12, data[[12]])
    
    sim_data_11 <- sim_dataset_fun(data[[11]])
    pop_pairs_11 <- pairs_creator(sim_data_11, data[[11]])
    fst_11 <- genet_dist_pairs(pop_pairs_11, method = "WC84")
    dch_11 <- genet_dist_pairs(pop_pairs_11, method = "Dch")
    jost_11 <- jostD_pairs(pop_pairs_11, data[[11]])
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    pop_pairs_10 <-pairs_creator(sim_data_10, data[[10]])
    fst_10 <- genet_dist_pairs(pop_pairs_10, method = "WC84")
    dch_10 <- genet_dist_pairs(pop_pairs_10, method = "Dch")
    jost_10 <- jostD_pairs(pop_pairs_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[09]])
    pop_pairs_09 <- pairs_creator(sim_data_09, data[[09]])
    fst_09 <- genet_dist_pairs(pop_pairs_09, method = "WC84")
    dch_09 <- genet_dist_pairs(pop_pairs_09, method = "Dch")
    jost_09 <- jostD_pairs(pop_pairs_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[08]])
    pop_pairs_08 <- pairs_creator(sim_data_08, data[[08]])
    fst_08 <- genet_dist_pairs(pop_pairs_08, method = "WC84")
    dch_08 <- genet_dist_pairs(pop_pairs_08, method = "Dch")
    jost_08 <- jostD_pairs(pop_pairs_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[07]])
    pop_pairs_07 <- pairs_creator(sim_data_07, data[[07]])
    fst_07 <- genet_dist_pairs(pop_pairs_07, method = "WC84")
    dch_07 <- genet_dist_pairs(pop_pairs_07, method = "Dch")
    jost_07 <- jostD_pairs(pop_pairs_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[06]])
    pop_pairs_06 <- pairs_creator(sim_data_06, data[[06]])
    fst_06 <- genet_dist_pairs(pop_pairs_06, method = "WC84")
    dch_06 <- genet_dist_pairs(pop_pairs_06, method = "Dch")
    jost_06 <- jostD_pairs(pop_pairs_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[05]])
    pop_pairs_05 <- pairs_creator(sim_data_05, data[[05]])
    fst_05 <- genet_dist_pairs(pop_pairs_05, method = "WC84")
    dch_05 <- genet_dist_pairs(pop_pairs_05, method = "Dch")
    jost_05 <- jostD_pairs(pop_pairs_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[04]])
    pop_pairs_04 <- pairs_creator(sim_data_04, data[[04]])
    fst_04 <- genet_dist_pairs(pop_pairs_04, method = "WC84")
    dch_04 <- genet_dist_pairs(pop_pairs_04, method = "Dch")
    jost_04 <- jostD_pairs(pop_pairs_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[03]])
    pop_pairs_03 <- pairs_creator(sim_data_03, data[[03]])
    fst_03 <- genet_dist_pairs(pop_pairs_03, method = "WC84")
    dch_03 <- genet_dist_pairs(pop_pairs_03, method = "Dch")
    jost_03 <- jostD_pairs(pop_pairs_03, data[[03]])
    
    sim_data_02 <- sim_dataset_fun(data[[02]])
    pop_pairs_02 <- pairs_creator(sim_data_02, data[[02]])
    fst_02 <- genet_dist_pairs(pop_pairs_02, method = "WC84")
    dch_02 <- genet_dist_pairs(pop_pairs_02, method = "Dch")
    jost_02 <- jostD_pairs(pop_pairs_02, data[[02]])
    
    
    
  }else if(id %in% c("Fagus_GR_Adult", "Fagus_GR_Regen", "Fagus_GR_Seed", 
                     "Fagus_SL_Adult", "Fagus_SL_Regen", "Fagus_SL_Seed")){
    
    sim_data_16 <- sim_dataset_fun(data[[16]])
    pop_pairs_16 <- pairs_creator(sim_data_16, data[[16]])
    fst_16 <- genet_dist_pairs(pop_pairs_16, method = "WC84")
    dch_16 <- genet_dist_pairs(pop_pairs_16, method = "Dch")
    jost_16 <- jostD_pairs(pop_pairs_16, data[[16]])
    
    sim_data_15 <- sim_dataset_fun(data[[15]])
    pop_pairs_15 <- pairs_creator(sim_data_15, data[[15]])
    fst_15 <- genet_dist_pairs(pop_pairs_15, method = "WC84")
    dch_15 <- genet_dist_pairs(pop_pairs_15, method = "Dch")
    jost_15 <- jostD_pairs(pop_pairs_15, data[[15]])
    
    sim_data_14 <- sim_dataset_fun(data[[14]])
    pop_pairs_14 <- pairs_creator(sim_data_14, data[[14]])
    fst_14 <- genet_dist_pairs(pop_pairs_14, method = "WC84")
    dch_14 <- genet_dist_pairs(pop_pairs_14, method = "Dch")
    jost_14 <- jostD_pairs(pop_pairs_14, data[[14]])
    
    sim_data_13 <- sim_dataset_fun(data[[13]])
    pop_pairs_13 <- pairs_creator(sim_data_13, data[[13]])
    fst_13 <- genet_dist_pairs(pop_pairs_13, method = "WC84")
    dch_13 <- genet_dist_pairs(pop_pairs_13, method = "Dch")
    jost_13 <- jostD_pairs(pop_pairs_13, data[[13]])
    
    sim_data_12 <- sim_dataset_fun(data[[12]])
    pop_pairs_12 <- pairs_creator(sim_data_12, data[[12]])
    fst_12 <- genet_dist_pairs(pop_pairs_12, method = "WC84")
    dch_12 <- genet_dist_pairs(pop_pairs_12, method = "Dch")
    jost_12 <- jostD_pairs(pop_pairs_12, data[[12]])
    
    sim_data_11 <- sim_dataset_fun(data[[11]])
    pop_pairs_11 <- pairs_creator(sim_data_11, data[[11]])
    fst_11 <- genet_dist_pairs(pop_pairs_11, method = "WC84")
    dch_11 <- genet_dist_pairs(pop_pairs_11, method = "Dch")
    jost_11 <- jostD_pairs(pop_pairs_11, data[[11]])
    
    sim_data_10 <- sim_dataset_fun(data[[10]])
    pop_pairs_10 <-pairs_creator(sim_data_10, data[[10]])
    fst_10 <- genet_dist_pairs(pop_pairs_10, method = "WC84")
    dch_10 <- genet_dist_pairs(pop_pairs_10, method = "Dch")
    jost_10 <- jostD_pairs(pop_pairs_10, data[[10]])
    
    sim_data_09 <- sim_dataset_fun(data[[09]])
    pop_pairs_09 <- pairs_creator(sim_data_09, data[[09]])
    fst_09 <- genet_dist_pairs(pop_pairs_09, method = "WC84")
    dch_09 <- genet_dist_pairs(pop_pairs_09, method = "Dch")
    jost_09 <- jostD_pairs(pop_pairs_09, data[[09]])
    
    sim_data_08 <- sim_dataset_fun(data[[08]])
    pop_pairs_08 <- pairs_creator(sim_data_08, data[[08]])
    fst_08 <- genet_dist_pairs(pop_pairs_08, method = "WC84")
    dch_08 <- genet_dist_pairs(pop_pairs_08, method = "Dch")
    jost_08 <- jostD_pairs(pop_pairs_08, data[[08]])
    
    sim_data_07 <- sim_dataset_fun(data[[07]])
    pop_pairs_07 <- pairs_creator(sim_data_07, data[[07]])
    fst_07 <- genet_dist_pairs(pop_pairs_07, method = "WC84")
    dch_07 <- genet_dist_pairs(pop_pairs_07, method = "Dch")
    jost_07 <- jostD_pairs(pop_pairs_07, data[[07]])
    
    sim_data_06 <- sim_dataset_fun(data[[06]])
    pop_pairs_06 <- pairs_creator(sim_data_06, data[[06]])
    fst_06 <- genet_dist_pairs(pop_pairs_06, method = "WC84")
    dch_06 <- genet_dist_pairs(pop_pairs_06, method = "Dch")
    jost_06 <- jostD_pairs(pop_pairs_06, data[[06]])
    
    sim_data_05 <- sim_dataset_fun(data[[05]])
    pop_pairs_05 <- pairs_creator(sim_data_05, data[[05]])
    fst_05 <- genet_dist_pairs(pop_pairs_05, method = "WC84")
    dch_05 <- genet_dist_pairs(pop_pairs_05, method = "Dch")
    jost_05 <- jostD_pairs(pop_pairs_05, data[[05]])
    
    sim_data_04 <- sim_dataset_fun(data[[04]])
    pop_pairs_04 <- pairs_creator(sim_data_04, data[[04]])
    fst_04 <- genet_dist_pairs(pop_pairs_04, method = "WC84")
    dch_04 <- genet_dist_pairs(pop_pairs_04, method = "Dch")
    jost_04 <- jostD_pairs(pop_pairs_04, data[[04]])
    
    sim_data_03 <- sim_dataset_fun(data[[03]])
    pop_pairs_03 <- pairs_creator(sim_data_03, data[[03]])
    fst_03 <- genet_dist_pairs(pop_pairs_03, method = "WC84")
    dch_03 <- genet_dist_pairs(pop_pairs_03, method = "Dch")
    jost_03 <- jostD_pairs(pop_pairs_03, data[[03]])
    
    sim_data_02 <- sim_dataset_fun(data[[02]])
    pop_pairs_02 <- pairs_creator(sim_data_02, data[[02]])
    fst_02 <- genet_dist_pairs(pop_pairs_02, method = "WC84")
    dch_02 <- genet_dist_pairs(pop_pairs_02, method = "Dch")
    jost_02 <- jostD_pairs(pop_pairs_02, data[[02]])
    
  }else{
    print("Unknown population - Cannot continue")
  }
}) 

# Plots ####  

pdf(paste(id, "distances_100_repl.pdf", sep = "_"), 
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

y_axis_fst <- seq(0, 0.999, 0.01)

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
