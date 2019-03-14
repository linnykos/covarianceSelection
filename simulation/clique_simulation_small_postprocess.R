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
hyp_tpr_list <- lapply(res, function(x){
  our_mat <- lapply(x, function(y){
    sapply(y$indices_list, function(z){
      length(which(idx %in% z))/length(idx_null)
    })
  })
  
  bonferroni_mat <- lapply(x, function(y){
    sapply(y$bonferroni_indices_list, function(z){
      length(which(idx %in% z))/length(idx_null)
    })
  })
  
  bh_mat <- lapply(x, function(y){
    sapply(y$bh_indices_list, function(z){
      length(which(idx %in% z))/length(idx_null)
    })
  })
  
  list(our_mat = our_mat, bonferroni_mat = bonferroni_mat, 
       bh_mat = bh_mat)
})

hyp_fpr_list <- lapply(res, function(x){
  our_mat <- lapply(x, function(y){
    sapply(y$indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })
  
  bonferroni_mat <- lapply(x, function(y){
    sapply(y$bonferroni_indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })
  
  bh_mat <- lapply(x, function(y){
    sapply(y$bh_indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })
  
  list(our_mat = our_mat, bonferroni_mat = bonferroni_mat, 
       bh_mat = bh_mat)
})

##

darkred <- rgb(205,40,54, max = 255)
             
png("../figures/individual_lvl1.png", height = 800, width = 2000, res = 300, units = "px")
trials <- length(hyp_fpr_list[[1]]$our_mat)
lvl <- 1
par(mfrow = c(1,3))
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Our method\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$our_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$our_mat[[i]],0), 
        col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (Bonferroni)\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$bonferroni_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$bonferroni_mat[[i]],0), 
        col = rgb(1, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (BH)\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$bh_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$bh_mat[[i]],0), 
        col = rgb(0.5, 0.5, 1, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
graphics.off()


png("../figures/individual_lvl2.png", height = 800, width = 2000, res = 300, units = "px")
trials <- length(hyp_fpr_list[[1]]$our_mat)
lvl <- 2
par(mfrow = c(1,3))
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Our method\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$our_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$our_mat[[i]],0), 
        col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (Bonferroni)\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$bonferroni_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$bonferroni_mat[[i]],0), 
        col = rgb(1, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (BH)\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))
for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$bh_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$bh_mat[[i]],0), 
        col = rgb(0.5, 0.5, 1, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
graphics.off()

#######################################

# now do the same but on the partition level

partition_tpr_list <- lapply(res, function(x){
  our_mat <- lapply(x, function(y){
    sapply(y$partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]
    })
  })
  
  bonferroni_mat <- lapply(x, function(y){
    sapply(y$bonferroni_partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]
    })
  })
  
  bh_mat <- lapply(x, function(y){
    sapply(y$bh_partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]
    })
  })
  
  list(our_mat = our_mat, bonferroni_mat = bonferroni_mat, 
       bh_mat = bh_mat)
})

partition_fpr_list <- lapply(res, function(x){
  our_mat <- lapply(x, function(y){
    sapply(y$partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(!z %in% 1:paramMat[1,1]))/(paramMat[1,2]+paramMat[1,3])
    })
  })
  
  bonferroni_mat <- lapply(x, function(y){
    sapply(y$bonferroni_partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(!z %in% 1:paramMat[1,1]))/(paramMat[1,2]+paramMat[1,3])
    })
  })
  
  bh_mat <- lapply(x, function(y){
    sapply(y$bh_partition_list, function(z){
      if(all(is.na(z))) return(0)
      length(which(!z %in% 1:paramMat[1,1]))/(paramMat[1,2]+paramMat[1,3])
    })
  })
  
  list(our_mat = our_mat, bonferroni_mat = bonferroni_mat, 
       bh_mat = bh_mat)
})


##

darkred <- rgb(205,40,54, max = 255)

png("../figures/partition_lvl1.png", height = 800, width = 2000, res = 300, units = "px")
trials <- length(partition_fpr_list[[1]]$our_mat)
lvl <- 1
par(mfrow = c(1,3))
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Our method\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$our_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$our_mat[[i]],0), 
        col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (Bonferroni)\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$bonferroni_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$bonferroni_mat[[i]],0), 
        col = rgb(1, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (BH)\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$bh_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$bh_mat[[i]],0), 
        col = rgb(0.5, 0.5, 1, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
graphics.off()




png("../figures/partition_lvl2.png", height = 800, width = 2000, res = 300, units = "px")
trials <- length(partition_fpr_list[[1]]$our_mat)
lvl <- 2
par(mfrow = c(1,3))
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Our method\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$our_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$our_mat[[i]],0), 
        col = rgb(0.5, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (Bonferroni)\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$bonferroni_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$bonferroni_mat[[i]],0), 
        col = rgb(1, 0.5, 0.5, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Naive method (BH)\nKappa = ", paramMat[lvl,"kappa"], ": Partitions"))
for(i in 1:trials){
  lines(c(1,partition_fpr_list[[lvl]]$bh_mat[[i]],0), c(1,partition_tpr_list[[lvl]]$bh_mat[[i]],0), 
        col = rgb(0.5, 0.5, 1, 0.5), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)
graphics.off()


