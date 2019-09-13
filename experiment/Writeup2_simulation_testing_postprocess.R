rm(list=ls())
load("../experiment/test_tmp.RData")
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 1000, c(0, 0.5, 1), seq(0, 1, length.out = 11)))
colnames(paramMat) <- c("num_group1", "num_group2", "num_group3", "n", "d",
                        "percentage", "alpha")

######

combn_null <- cbind(combn(paramMat[1,1],2),
                    (combn(paramMat[1,2],2)+paramMat[1,1]),
                    (combn(paramMat[1,3],2)+sum(paramMat[1,1:2])))
num_partition <- sum(paramMat[1,1:3])
idx_null <- combn_null[1,]+num_partition*combn_null[2,]
combn_mat <- combn(num_partition,2)
idx_all <- combn_mat[1,]+num_partition*combn_mat[2,]
idx <- which(idx_all %in% idx_null)

###########

# individual hypothesis
hyp_tpr_list <- lapply(res, function(x){
  sapply(x, function(y){
    length(which(idx %in% y$res))/length(idx_null)
  })
})

hyp_fpr_list <- lapply(res, function(x){
  sapply(x, function(y){
    length(which(!y$res %in% idx))/(length(idx_all) - length(idx_null))
  })
})

# reformat
lis_tpr <- vector("list", 3)
for(i in 1:3){
  mat <- matrix(NA, nrow = 11, ncol = 10)
  for(j in 1:11){
    mat[j,] <- hyp_tpr_list[[(j-1)*3+i]]
  }
  
  lis_tpr[[i]] <- mat
}

lis_fpr <- vector("list", 3)
for(i in 1:3){
  mat <- matrix(NA, nrow = 11, ncol = 10)
  for(j in 1:11){
    mat[j,] <- hyp_fpr_list[[(j-1)*3+i]]
  }
  
  lis_fpr[[i]] <- mat
}

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:10){
  yvec <- c(1,lis_tpr[[1]][,i],0)
  xvec <- c(1,lis_fpr[[1]][,i],0)
  lines(xvec, yvec)
}

#################

# let's try out subset selection method 
clique_func_tmp <- function(vec, len = 25){
  vec <- vec$res
  if(length(vec) == 0) return(NA)
  g <- igraph::graph.empty(n = len, directed = F)
  g <- igraph::add_edges(g, edges = combn(len, 2)[,vec])
  lis <- clique_selection(g)
  
  if(length(lis) > 1){
    intersect_len <- sapply(lis, function(x){length(which(1:15 %in% x))})
    lis[[which.max(intersect_len)]]
  } else{
    lis[[1]]
  }
}

res_partitions <- lapply(res, function(x){
  lapply(x, clique_func_tmp)
})

######################

partition_tpr_list <- lapply(res_partitions, function(x){
  sapply(x, function(y){
    length(which(1:paramMat[1,1] %in% y))/paramMat[1,1]
  })
})

partition_fpr_list <- lapply(res_partitions, function(x){
  sapply(x, function(y){
    length(which(!y %in% 1:paramMat[1,1]))/(paramMat[1,2]+paramMat[1,3])
  })
})

# reformat
lis_partition_tpr <- vector("list", 3)
for(i in 1:3){
  mat <- matrix(NA, nrow = 11, ncol = 10)
  for(j in 1:11){
    mat[j,] <- partition_tpr_list[[(j-1)*3+i]]
  }
  
  lis_partition_tpr[[i]] <- mat
}

lis_partition_fpr <- vector("list", 3)
for(i in 1:3){
  mat <- matrix(NA, nrow = 11, ncol = 10)
  for(j in 1:11){
    mat[j,] <- partition_fpr_list[[(j-1)*3+i]]
  }
  
  lis_partition_fpr[[i]] <- mat
}

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:10){
  yvec <- c(1,lis_partition_tpr[[3]][,i],0)
  xvec <- c(1,lis_partition_fpr[[3]][,i],0)
  lines(xvec, yvec)
}
