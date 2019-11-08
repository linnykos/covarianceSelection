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

# remove res that errored
for(i in 1:length(res)){
  rm_idx <- which(sapply(res[[i]], length) == 1)
  if(length(rm_idx) > 0) res[[i]] <- res[[i]][-rm_idx]
}

# assign colors
colfunc <- colorRampPalette(c(rgb(205,40,54, max = 255), rgb(149,219,144, max = 255)))
col_vec <- colfunc(4)
lwd_vec <- c(5.5, 5, 4.5, 4)

###########

# individual hypothesis
hyp_tpr_list <- lapply(res, function(x){ #loop over simulation settings
  sapply(x, function(y){ #loop over trials
    c(1, sapply(1:21, function(i){ #loop over alphas
      length(which(idx %in% y$indices_list[[i]]))/length(idx_null)
    }), 0)
  })
})

hyp_fpr_list <- lapply(res, function(x){
  sapply(x, function(y){
    c(1, sapply(1:21, function(i){
      length(which(!y$indices_list[[i]] %in% idx))/(length(idx_all) - length(idx_null))
    }), 0)
  })
})

##############

# individual hypothesis
hyp_tpr_list2 <- lapply(res, function(x){ #loop over simulation settings
  sapply(x, function(y){ #loop over trials
    c(1, sapply(1:21, function(i){ #loop over alphas
      length(which(idx %in% y$bonferroni_indices_list[[i]]))/length(idx_null)
    }), 0)
  })
})

hyp_fpr_list2 <- lapply(res, function(x){
  sapply(x, function(y){
    c(1, sapply(1:21, function(i){
      length(which(!y$bonferroni_indices_list[[i]] %in% idx))/(length(idx_all) - length(idx_null))
    }), 0)
  })
})


# plot stepdown pvalues
png("../figures/figure_5.png", height = 1400, width = 2400, res = 300, units ="px")
par(mfrow = c(1,2), mar = c(5,4,4,1))

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T, xlab = "False positive rate", ylab = "True positive rate",
     main = "Individual hypotheses using\nnaive family−wise correction")
for(k in 1:4){
  roc_our <- roc_region(hyp_tpr_list2[[k]], hyp_fpr_list2[[k]])
  lines(c(0, seq(0, 1, length.out = 21), 1), c(0, roc_our$median_tpr, 1), col = col_vec[k], lwd = lwd_vec[k])
}
legend("bottomright", rev(c("0% level", "30% level", "60% level", "100% level")),
       bty="n", fill=rev(col_vec))

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T, xlab = "False positive rate", ylab = "True positive rate",
     main = "Individual hypotheses using\nour Stepdown method")
for(k in 1:4){
  roc_our <- roc_region(hyp_tpr_list[[k]], hyp_fpr_list[[k]])
  lines(c(0, seq(0, 1, length.out = 21), 1), c(0, roc_our$median_tpr, 1), col = col_vec[k], lwd = lwd_vec[k])
}
legend("bottomright", rev(c("0% level", "30% level", "60% level", "100% level")),
       bty="n", fill=rev(col_vec))
graphics.off()