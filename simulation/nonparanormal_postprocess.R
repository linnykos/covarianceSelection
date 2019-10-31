rm(list = ls())
load("../results/nonparanormal.RData")

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
hyp_tpr_list <- lapply(res, function(x){ #loop over simulation settings
  sapply(x, function(y){ #loop over trials
    sapply(1:21, function(i){ #loop over alphas
      length(which(idx %in% y$indices_list[[i]]))/length(idx_null)
    })
  })
})

hyp_fpr_list <- lapply(res, function(x){
  sapply(x, function(y){
    sapply(1:21, function(i){
      length(which(!y$indices_list[[i]] %in% idx))/(length(idx_all) - length(idx_null))
    })
  })
})

# plot stepdown pvalues
k <- 1
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:5){
  lines(c(1,hyp_fpr_list[[k]][,i],0), c(1,hyp_tpr_list[[k]][,i],0))
}

##############

# individual hypothesis
hyp_tpr_list <- lapply(res, function(x){ #loop over simulation settings
  sapply(x, function(y){ #loop over trials
    sapply(1:21, function(i){ #loop over alphas
      length(which(idx %in% y$bh_indices_list[[i]]))/length(idx_null)
    })
  })
})

hyp_fpr_list <- lapply(res, function(x){
  sapply(x, function(y){
    sapply(1:21, function(i){
      length(which(!y$bh_indices_list[[i]] %in% idx))/(length(idx_all) - length(idx_null))
    })
  })
})

# plot stepdown pvalues
k <- 3
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:5){
  lines(c(1,hyp_fpr_list[[k]][,i],0), c(1,hyp_tpr_list[[k]][,i],0))
}

