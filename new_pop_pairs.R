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




sim_data_01 <- sim_dataset_fun(data[[01]])
sim_dataset <- sim_data_01
empirical <- data[[01]]
########################
# The following code does not throw warnings!


pop_name <- list()
temp <- list()
for(i in samp_size[-length(samp_size)]){
  for(j in 1:replic_num){
    pop_name[[i]][[j]] <- list(paste("pop", i, j, sep = "_"))
    pop(sim_dataset[[i]][[j]]) <- rep(pop_name[[i]][[j]], nrow(sim_dataset[[i]][[j]]@tab))
  }
}

pop_name[[samp_size[length(samp_size)]]] <- paste("pop", samp_size[length(samp_size)], sep = "_")

pop(sim_dataset[[samp_size[length(samp_size)]]]) <- 
  rep(pop_name[[samp_size[length(samp_size)]]], nrow(sim_dataset[[samp_size[length(samp_size)]]]@tab) )


Allpops <- repool(unlist(sim_dataset)) # Create a total pool
Allpops <- repool(empirical, Allpops)

temp <- list()
for(i in samp_size[-length(samp_size)]){
  for(j in 1:replic_num){
    temp[[i]][[j]] <- list(Allpops[pop = c(paste(pop_name[[i]][[j]]), "GR_Seed")])
  }
}


temp[[samp_size[length(samp_size)]]] <- 
  Allpops[pop = c(paste(pop_name[[samp_size[length(samp_size)]]]), "GR_Seed")]
