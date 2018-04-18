rm(list=ls())
load("../results/simulation_RoC.RData")
load("../results/simulation_RoC_part2.RData")

par(mfrow = c(1,3))

plot(NA, xlim = c(1, 21), ylim = c(0, 25))
for(i in 1:20){
  for(j in 4:4){
    tmp <- sapply(res_indices_list_spectral[[i]][[j]], length)
    lines(tmp, col = rgb(0.2, 0.2, 0.2, 0.2))
  }
}

plot(NA, xlim = c(1, 21), ylim = c(0, 25))
for(i in 1:20){
  for(j in 4:4){
    tmp <- sapply(res_indices_list_subgraph_avg[[i]][[j]], length)
    lines(tmp, col = rgb(0.2, 0.2, 0.2, 0.2))
  }
}

plot(NA, xlim = c(1, 21), ylim = c(0, 25))
for(i in 1:20){
  for(j in 4:4){
    tmp <- sapply(res_partition_len_list[[i]][[j]], length)
    lines(tmp, col = rgb(0.2, 0.2, 0.2, 0.2))
  }
}
