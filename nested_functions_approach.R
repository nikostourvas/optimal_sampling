# Original function
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

# Nested function
sims <- function(data){
  res <- list()
  for(j in 1:length(data)){
    res[[j]] <- sim_dataset_fun(data[[j]])
  }
  return(res)  
}  


x <- sims(data)
