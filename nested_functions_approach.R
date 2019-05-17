
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


# Hexp
extract_het_fun <- function(x){
  x$Hobs
}


hexp_rap <- rapply(res_rap, extract_het_fun, how = "replace")


# ---

library(rlang)
test5 <- squash(res_rap)
test6 <- melt(
  lapply(test5, extract_het_fun)
)

repl2 <- rep(1:10, each = 2, times = 10)
repl3 <- c(repl2, 1)
repl4 <- rep(repl3, nLoc(obj))

# ---

# Original function
Het_fun <- function(sim_dataset, results, heterozygosity){
  # Ho for each generated dataset 
  Het <- list()
  for(i in samp_size[-length(samp_size)]){
    for(j in 1:replic_num){ # number of replications
      Het[[i]][[j]] <- results[[i]][[j]][[heterozygosity]]
    }
  }
  
  Het[[samp_size[length(samp_size)]]] <- 
    results[[samp_size[length(samp_size)]]][[heterozygosity]]
  
  # Mean Ho values for each generated dataset 
  Het_means <- list()
  for(i in samp_size[-length(samp_size)]){
    Het_means[[i]] <- lapply(Het[[i]], mean)
  }
  
  Het_means[[samp_size[length(samp_size)]]] <- 
    mean(Het[[samp_size[length(samp_size)]]])
  
  # Produce a data frame to be plotted by ggplot2
  Het_means_df <- melt(Het_means)
  colnames(Het_means_df) <- c("value", "replic", "samp_size")
  Het_means_df$samp_size <- 
    factor(Het_means_df$samp_size, levels = unique(
      as.character(Het_means_df$samp_size)))
  
  Het_means_df$marker_num <- 
    if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
    paste(length(nAll(sim_dataset[[1]][[1]])), 
          "markers", sep = " ")
  }else {
    paste(length(nAll(sim_dataset[[1]][[1]])), 
          "marker", sep = " ")
  }
  
  return(Het_means_df)
}  


# Nested function
het_fun <- function(sim_dataset){
  het <- list()
  for(j in 1:length(sim_dataset)){
    het[[j]] <- Het_fun(sim_dataset[[j]], res_rap, "Hobs")
  }
  return(het)
}

het <- het_fun(sims)

