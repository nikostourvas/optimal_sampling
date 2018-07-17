empirical <- read.genalexcel(
  "LGM_DE_SI_GR_final.xlsx",   # name of excel file
  sheet = "Abies_emp",             # name of sheet where the genotypes reside
  genclone = F)

empirical <- missingno(empirical, type = "mean")

pop <- "emp_GR_Adult" # select pop to analyze (acceptable names: "DE_Adult", ")

# Seperate dataset by pop & select pop to analyze
emp_list <- seppop(empirical) # separate pops
empirical <- emp_list[[pop]]

data <- list()
for(i in 1:length(loci)){
  data[[length(loci)+1-i]] <- empirical[, loc = loci[i:length(loci)]]
}

#####################################################################################
Dch_pairs <- function(sim_dataset, empirical){
  
  # Rename pop of empirical dataset, so that it can be distinguished from the replicate.
  data_emp <- empirical
  pop(data_emp) <- rep("emp", nrow(empirical@tab))
  
      
  # Create list with one replicate & the empirical dataset
  pop_pairs_list <- list()
  for(i in samp_size){
    for(j in 1:replic_num){
      pop_pairs_list[[i]][[j]] <- repool(sim_dataset[[i]][[j]], data_emp, list = TRUE)
    }
  }
  
  # Create genind objects with one replicate & the empirical dataset
  pop_pairs <- list()
  for(i in samp_size){
    pop_pairs[[i]] <- lapply(pop_pairs_list[[i]], repool)
  }
  
  # Calculate distance
  Cav_Sf <- list()
  for(i in samp_size){
    Cav_Sf[[i]] <- as.data.frame(sapply(pop_pairs[[i]], genet.dist, method = "Dch"))
  }
  
  for(i in samp_size){
    Cav_Sf[[i]]$samp_size <- paste(i)
    colnames(Cav_Sf[[i]]) <- c("Dch", "samp_size")
  }
  
  # Create a single data.frame with all the distances
  Cav_Sf_df <- as.data.frame(bind_rows(Cav_Sf))
  Cav_Sf_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
                              paste(length(nAll(sim_dataset[[1]][[1]])), 
                                    "markers", sep = " ")
                            }else {
                              paste(length(nAll(sim_dataset[[1]][[1]])), 
                                    "marker", sep = " ")
                            }

  return(Cav_Sf_df)
}

test <- Dch_pairs(sim_data_10, data[[10]])
############################################################################
# Rename pop of empirical dataset, so that it can be distinguished from the replicate.
data_emp <- data
pop(data_emp[[2]]) <- rep("emp_GR_Adult", nrow(data[[2]]@tab))

test <- repool(sim_data_02[["50"]][[1]], data_emp[[2]])
dch <- genet.dist(test)
