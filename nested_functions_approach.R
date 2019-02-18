
# Simulations -------------------------------------------------------------


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
sim_fun <- function(data){
  res <- list()
  for(j in 1:length(data)){
    res[[j]] <- sim_dataset_fun(data[[j]])
  }
  return(res)  
}  


sims <- sim_fun(data)


# Hobs, Hexp --------------------------------------------------------------

# Original function
results_fun <- function(sim_dataset){
  results <- list()
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){ # number of replications
      results[[i]][[j]] <- summary(sim_dataset[[i]][[j]])
    }
  }
  
  results[[samp_size[length(samp_size)]]] <- 
    summary(sim_dataset[[samp_size[length(samp_size)]]])
  
  return(results)
}

# Nested function
res_fun <- function(sim_dataset){
  res <- list()
  for(j in 1:length(sim_dataset)){
    res[[j]] <- results_fun(sim_dataset[[j]])
  }
  return(res)
}

res <- res_fun(sims)

# rapply way
# not faster, but way more elegant

system.time({
res_rap <- rapply(sims, summary, how = "replace")
})

system.time({
  res <- res_fun(sims)
})
