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

############################

tpr_mat <- rbind(1, do.call(cbind, hyp_tpr_list[[1]]$our_mat), 0)
fpr_mat <- rbind(1, do.call(cbind, hyp_fpr_list[[1]]$our_mat), 0)

fpr_seq <- seq(0, 1, length.out = 51)
res <- roc_region(tpr_mat, fpr_mat, fpr_seq = fpr_seq, quant = c(0,1))

#############################

lvl <- 1
trials <- 20
darkred <- rgb(205,40,54, max = 255)
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate", 
     asp = T, main = paste0("Our method\nKappa = ", paramMat[lvl,"kappa"], ": Pairs"))

polygon(x = c(fpr_seq, rev(fpr_seq)), y = c(res$lower_tpr, rev( res$upper_tpr)),
        col = rgb(0,0,1,0.2), lwd = 3,
        border = rgb(0,0,1))
lines(fpr_seq, res$median_tpr, lwd = 5, col = "blue")

for(i in 1:trials){
  lines(c(1,hyp_fpr_list[[lvl]]$our_mat[[i]],0), c(1,hyp_tpr_list[[lvl]]$our_mat[[i]],0), 
        col = rgb(0.5, 0.5, 0.5, 0.2), lwd = 2)
}
lines(c(0,1), c(0,1), col = darkred, lty = 2, lwd = 4)

# plot the region
# lines(fpr_seq, res$lower_tpr, lwd = 5, lty = 2, col = "blue")
# lines(fpr_seq, res$upper_tpr, lwd = 5, lty = 2, col = "blue")
# lines(fpr_seq, res$median_tpr, lwd = 5, col = "blue")

