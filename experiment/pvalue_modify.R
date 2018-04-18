rm(list=ls())
load("../results/simulation_RoC.RData")

alpha_vec <- seq(0, 1, length.out = 21)

for(i in 1:4){
  for(j in 1:50){
    res[[i]][[j]]$naive_indices_list <- lapply(alpha_vec, function(alpha){
      which(res[[i]][[j]]$naive_pval_vec>= alpha)
    })
  }
}

save.image("../results/simulation_RoC.RData")
