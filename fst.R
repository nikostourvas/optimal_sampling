empirical <- read.genalexcel(
  "LGM_DE_SI_GR_final.xlsx",   # name of excel file
  sheet = "Abies_emp",             # name of sheet where the genotypes reside
  genclone = F)

seppop(empirical)
emp_list <- seppop(empirical)
DE_Adult <- emp_list[[1]]
GR_Adult <- emp_list[[7]]
GR_Regen <- emp_list[[8]]


pop_pairs <- list()
for(i in samp_size){
  for(j in 1:replic_num){
    pop_pairs[[i]][[j]] <- repool(sim_data_11[[i]][[j]], data[[11]], list = TRUE)
  }
}

pop_pairs_200 <- lapply(pop_pairs[["200"]], repool)
Cav_Sf_200 <- as.data.frame(sapply(pop_pairs_200, genet.dist, method = "Dch"))

pop_pairs_test <- list()
for(i in samp_size){
  pop_pairs_test[[i]] <- lapply(pop_pairs[[i]], repool)
}

Cav_Sf <- list()
for(i in samp_size){
  Cav_Sf[[i]] <- as.data.frame(sapply(pop_pairs_test[[i]], genet.dist, method = "Dch"),
                               colnames = i)
}








Cav_Sf <- list()
for(i in samp_size){
  for(j in 1:replic_num){
Cav_Sf[[i]][[j]] <- genet.dist(pop_pairs[[i]][[j]], method = "Dch")
    }
}


DS_250 <- as.data.frame(0)

pop_pairs_225 <- list()
for(j in 1:replic_num){
  pop_pairs_225[[j]] <- repool(data[[11]], sim_data_11[["225"]][[j]])
}
DS_225 <- as.data.frame(sapply(pop_pairs_225, genet.dist, method = "Dch"))
DS_225$samp_size <- "225"

pop_pairs_200 <- list()
  for(j in 1:replic_num){
    pop_pairs_200[[j]] <- repool(data[[11]], sim_data_11[["200"]][[j]])
  }
DS_200 <- as.data.frame(sapply(pop_pairs_200, genet.dist, method = "Dch"))
DS_200$samp_size <- "200"

pop_pairs_175 <- list()
for(j in 1:replic_num){
  pop_pairs_175[[j]] <- repool(data[[11]], sim_data_11[["175"]][[j]])
}
DS_175 <- as.data.frame(sapply(pop_pairs_175, genet.dist, method = "Dch"))
DS_175$samp_size <- "175"







DS <- list()
for(j in 1:replic_num){
DS[[j]] <- genet.dist(pop_pairs_250[[j]], method ="Dch")
}




DS <- lapply(pop_pairs_250, genet.dist)

pegas_DS <- list()
for(j in 1:replic_num){
pegas_DS[[j]] <- Fst(as.loci(pop_pairs_250[[j]]))
}


pop_pairs <- lapply(data[[11]], sim_data_11, repool)

pop_pairs <- repool(data[[11]], sim_data_11[[1]][[1]])





                    