# Convert genind files to fstat file format
# to calculate Ne in Ne Estimator

# Install (if needed) & load required libraries
if (!require("adegenet")) install.packages("adegenet")
library(adegenet)
if (!require("popprxl")) install.packages("popprxl")
library(popprxl)
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

# Setup -------------------------------------------------------------------


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


# hierfstat ---------------------------------------------------------------

library(hierfstat)

hier_list <- rapply(sims, genind2hierfstat, how = "replace")


for(h in 2:11){
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){ # number of replications
      write.fstat(hier_list[[h]][[i]][[j]],
                  fname = paste(h, i, j, "obj.dat", sep = "_"))
    }
  }
}

for(h in 2:11){
write.fstat(hier_list[[h]][[samp_size[length(samp_size)]]],
            fname = paste(h, samp_size[length(samp_size)], 
                                       "obj.dat", sep="_"))
}


# write.fstat(obj_hier,
#             fname = paste(i, "test.dat", sep = "_") )
