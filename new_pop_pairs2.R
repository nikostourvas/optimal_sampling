# first create pop_pairs_all

pairs_initial <- function(sim_dataset){
  # Rename pop of empirical dataset, so that it can be distinguished from the replicate.
  pop(data[[length(loci)]]) <- rep("emp", nrow(data[[length(loci)]]@tab))
  
  # Create list with one replicate & the empirical dataset
  pop_pairs_list <- list()
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){
      pop_pairs_list[[i]][[j]] <- repool(sim_dataset[[i]][[j]], data[[length(loci)]], list = TRUE)
    }
  }
  
  pop_pairs_list[[samp_size[length(samp_size)]]] <- 
    repool(sim_dataset[[samp_size[length(samp_size)]]], data[[length(loci)]], list = FALSE)
  
  
  # Create genind objects with one replicate & the empirical dataset
  pop_pairs_all <- list()
  for(i in samp_size[-length(samp_size)]){
    pop_pairs_all[[i]] <- lapply(pop_pairs_list[[i]], repool)
  }
  
  pop_pairs_all[[samp_size[length(samp_size)]]] <- pop_pairs_list[[samp_size[length(samp_size)]]]
  
  
  return(pop_pairs_all)
}
#################################################################
pairs_creator <- function(marker_num){

  x <- length(loci) - marker_num + 1
  pop_pairs <- list()
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){  
      pop_pairs[[i]][[j]] <- list(pop_pairs_all[[i]][[j]][, loc = loci[x:length(loci)]])
    }
  }
  
  for(i in samp_size[-length(samp_size)]){
    pop_pairs[[i]] <- unlist(pop_pairs[[i]])
  }
  
  pop_pairs[[samp_size[length(samp_size)]]] <- 
    pop_pairs_all[[samp_size[length(samp_size)]]][, loc = loci[x:length(loci)]]
  
  return(pop_pairs)
}

##############################################
for(i in samp_size[-length(samp_size)]){
  pop_pairs[[i]] <- unlist(pop_pairs[[i]])
}
