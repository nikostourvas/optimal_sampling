input <- data[[10]]
sim_dataset <- sim_data
empirical <- data[[10]]
method <- "WC84"

# functions

# sim_dataset_fun without creating multiple iterations of the empirical data set
sim_dataset_fun <- function(input){
  sim_data <- list()
  for(i in samp_size[-length(samp_size)]){ 
    sim_data[[i]] <-
      replicate (replic_num, input[sample(1:nrow(input$tab), 
                                          i, replace = F)])
  }
  
  sim_data[[samp_size[length(samp_size)]]] <- input
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



genet_dist_pairs <- function(sim_dataset, empirical, method){

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
  distance_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){ 
    paste(length(nAll(sim_dataset[[1]][[1]])), 
          "markers", sep = " ")
  }else {
    paste(length(nAll(sim_dataset[[1]][[1]])), 
          "marker", sep = " ")
  }
  
  # Replace negative values with zero
  distance_df[, method][distance_df[, method] < 0] <- 0
  
  
  return(distance_df)
}




jostD_pairs <- function(sim_dataset, empirical){
  
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
  distance_df$marker_num <- if(length(nAll(sim_dataset[[1]][[1]])) > 1){
    paste(length(nAll(sim_dataset[[1]][[1]])),
          "markers", sep = " ")
  }else {
    paste(length(nAll(sim_dataset[[1]][[1]])),
          "marker", sep = " ")
  }
  
  # Replace negative values with zero
  distance_df[, "D_Jost"][distance_df[, "D_Jost"] < 0] <- 0
  
  return(distance_df)
}