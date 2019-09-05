rm(list=ls())
load("../results/clique_simulation_small.RData")

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
order_idx <- c("indices_list", "bonferroni_indices_list", "bh_indices_list")
hyp_tpr_list <- lapply(res, function(x){
  lis <- vector("list", length(order_idx))
  for(i in 1:length(order_idx)){
    lis[[i]] <- sapply(x, function(y){
      j <- which(names(y) == order_idx[i])
      c(1, sapply(y[[j]], function(z){
        length(which(idx %in% z))/length(idx_null)
      }), 0)
    })
  }
  
  names(lis) = c("our_mat", "bonferroni_mat", "bh_mat")
  
  lis
})

hyp_fpr_list <- lapply(res, function(x){
  lis <- vector("list", length(order_idx))
  for(i in 1:length(order_idx)){
    lis[[i]] <- sapply(x, function(y){
      j <- which(names(y) == order_idx[i])
      c(1, sapply(y[[j]], function(z){
        length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
      }), 0)
    })
  }
  
  names(lis) = c("our_mat", "bonferroni_mat", "bh_mat")
  
  lis
})

##

darkred <- rgb(205,40,54, max = 255)
             
png("../figures/individual_lvl1.png", height = 800, width = 2000, res = 300, units = "px")
trials <- ncol(hyp_fpr_list[[1]][[1]])
lvl <- 1
name_vec <- c("Our method", "Naive method (Bonferroni)", "Naive method (BH)")
par(mfrow = c(1,length(name_vec)))
for(j in 1:length(name_vec)){
  fpr_seq <- seq(0, 1, length.out = 51)
  res_roc <- roc_region(hyp_tpr_list[[lvl]][[j]], hyp_fpr_list[[lvl]][[j]], fpr_seq = fpr_seq, quant = c(0,1))
  
  plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
       asp = T, main = paste0(name_vec[j], "\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
  
  polygon(x = c(fpr_seq, rev(fpr_seq)), y = c(res_roc$lower_tpr, rev(res_roc$upper_tpr)),
          col = rgb(0,0,1,0.2), lwd = 2,
          border = rgb(0,0,1))
  
  
  for(i in 1:trials){
    lines(c(1,hyp_fpr_list[[lvl]][[j]][,i],0), c(1,hyp_tpr_list[[lvl]][[j]][,i],0), 
          col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 1)
  }
  lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
  
  lines(c(0,fpr_seq,1), c(0,res_roc$median_tpr,1), lwd = 5, col = "blue")
  
}
graphics.off()


png("../figures/individual_lvl2.png", height = 800, width = 2000, res = 300, units = "px")
trials <- ncol(hyp_fpr_list[[1]][[1]])
lvl <- 2
name_vec <- c("Our method", "Naive method (Bonferroni)", "Naive method (BH)")
par(mfrow = c(1,length(name_vec)))
for(j in 1:length(name_vec)){
  fpr_seq <- seq(0, 1, length.out = 51)
  res_roc <- roc_region(hyp_tpr_list[[lvl]][[j]], hyp_fpr_list[[lvl]][[j]], fpr_seq = fpr_seq, quant = c(0,1))
  
  plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
       asp = T, main = paste0(name_vec[j], "\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
  
  polygon(x = c(fpr_seq, rev(fpr_seq)), y = c(res_roc$lower_tpr, rev(res_roc$upper_tpr)),
          col = rgb(0,0,1,0.2), lwd = 2,
          border = rgb(0,0,1))
  
  
  for(i in 1:trials){
    lines(c(1,hyp_fpr_list[[lvl]][[j]][,i],0), c(1,hyp_tpr_list[[lvl]][[j]][,i],0), 
          col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 1)
  }
  lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
  
  lines(c(0,fpr_seq,1), c(0,res_roc$median_tpr,1), lwd = 5, col = "blue")
  
}
graphics.off()

#######################################

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

##

darkred <- rgb(205,40,54, max = 255)

png("../figures/partition_lvl1.png", height = 800, width = 2000, res = 300, units = "px")
trials <- ncol(partition_fpr_list[[1]][[1]])
lvl <- 1
name_vec <- c("Our method", "Naive method (Bonferroni)", "Naive method (BH)")
par(mfrow = c(1,length(name_vec)))
for(j in 1:length(name_vec)){
  fpr_seq <- seq(0, 1, length.out = 51)
  res_roc <- roc_region(partition_tpr_list[[lvl]][[j]], partition_fpr_list[[lvl]][[j]], fpr_seq = fpr_seq, 
                        quant = c(0.1, 0.9))
  
  plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
       asp = T, main = paste0(name_vec[j], "\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
  
  polygon(x = c(fpr_seq, rev(fpr_seq)), y = c(res_roc$lower_tpr, rev(res_roc$upper_tpr)),
          col = rgb(0,0,1,0.2), lwd = 2,
          border = rgb(0,0,1))
  
  
  for(i in 1:trials){
    lines(c(1,partition_fpr_list[[lvl]][[j]][,i],0), c(1,partition_tpr_list[[lvl]][[j]][,i],0), 
          col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 1)
  }
  lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
  
  lines(c(0,fpr_seq,1), c(0,res_roc$median_tpr,1), lwd = 5, col = "blue")
  
}
graphics.off()


png("../figures/partition_lvl2.png", height = 800, width = 2000, res = 300, units = "px")
trials <- ncol(partition_fpr_list[[1]][[1]])
lvl <- 2
name_vec <- c("Our method", "Naive method (Bonferroni)", "Naive method (BH)")
par(mfrow = c(1,length(name_vec)))
for(j in 1:length(name_vec)){
  fpr_seq <- seq(0, 1, length.out = 51)
  res_roc <- roc_region(partition_tpr_list[[lvl]][[j]], partition_fpr_list[[lvl]][[j]], fpr_seq = fpr_seq, 
                        quant = c(0.1, 0.9))
  
  plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
       asp = T, main = paste0(name_vec[j], "\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
  
  polygon(x = c(fpr_seq, rev(fpr_seq)), y = c(res_roc$lower_tpr, rev(res_roc$upper_tpr)),
          col = rgb(0,0,1,0.2), lwd = 2,
          border = rgb(0,0,1))
  
  
  for(i in 1:trials){
    lines(c(1,partition_fpr_list[[lvl]][[j]][,i],0), c(1,partition_tpr_list[[lvl]][[j]][,i],0), 
          col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 1)
  }
  lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
  
  lines(c(0,fpr_seq,1), c(0,res_roc$median_tpr,1), lwd = 5, col = "blue")
  
}
graphics.off()