rm(list=ls())
load("../results/clique_simulation_small.RData")
# now do the same but on the partition level
order_idx <- c("partition_list", "bonferroni_partition_list", "bh_partition_list")
partition_tpr_list <- lapply(res, function(x){
  lis <- vector("list", length(order_idx))
  for(i in 1:length(order_idx)){
    lis[[i]] <- sapply(x, function(y){
      j <- which(names(y) == order_idx[i])
      c(1, sapply(y[[j]], function(z){
        if(all(is.na(z))) return(0)
        length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]
      }), 0)
    })
  }
  
  names(lis) = c("our_mat", "bonferroni_mat", "bh_mat")
  
  lis
})

partition_fpr_list <- lapply(res, function(x){
  lis <- vector("list", length(order_idx))
  for(i in 1:length(order_idx)){
    lis[[i]] <- sapply(x, function(y){
      j <- which(names(y) == order_idx[i])
      c(1, sapply(y[[j]], function(z){
        if(all(is.na(z))) return(0)
        length(which(!z %in% 1:paramMat[1,1]))/(paramMat[1,2]+paramMat[1,3])
      }), 0)
    })
  }
  
  names(lis) = c("our_mat", "bonferroni_mat", "bh_mat")
  
  lis
})