#plot the ROC curves

rm(list=ls())
load("../results/simulation_RoC.RData")

#####################

num_partition <- sum(paramMat[1,1:3])
stopifnot(length(unique(paramMat[,1])) == 1, length(unique(paramMat[,2])) == 1, length(unique(paramMat[,3])) == 1)

# for indices, compute TPR and FPR
combn_null <- cbind(combn(paramMat[1,1],2),
                    (combn(paramMat[1,2],2)+paramMat[1,1]),
                    (combn(paramMat[1,3],2)+sum(paramMat[1,1:2])))
idx_null <- combn_null[1,]+num_partition*combn_null[2,]
combn_mat <- combn(num_partition,2)
idx_all <- combn_mat[1,]+num_partition*combn_mat[2,]

# true postive rate (our method)
idx <- which(idx_all %in% idx_null)
hyp_tpr_list <- lapply(res, function(x){ #on each flip-rate
  mat <- sapply(x, function(y){ #trials
    sapply(y$indices_list, function(z){ #over alpha
      length(which(idx %in% z))/length(idx_null)
    })
  })

  apply(mat, 1, mean)
})

# true postive rate (naive method)
hyp_tpr_list_naive <- lapply(res, function(x){ #on each flip-rate
  mat <- sapply(x, function(y){ #trials
    sapply(y$naive_indices_list, function(z){ #over alpha
      length(which(idx %in% z))/length(idx_null)
    })
  })

  apply(mat, 1, mean)
})

# false positive rate (our method)
hyp_fpr_list <- lapply(res, function(x){
  mat <- sapply(x, function(y){
    sapply(y$indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })

  apply(mat, 1, mean)
})

# false positive rate (naive method)

hyp_fpr_list_naive <- lapply(res, function(x){
  mat <- sapply(x, function(y){
    sapply(y$naive_indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })

  apply(mat, 1, mean)
})



#######################

# create color ramp
colfunc <- colorRampPalette(c(rgb(205,40,54, max = 255), rgb(149,219,144, max = 255)))
col_vec <- colfunc(4)
lwd_vec <- c(5.5, 5, 4.5, 4)

#######################

alpha_vec <- seq(0, 1, length.out = 21)
alpha_vec <- round(alpha_vec*100)/100

pdf("../figures/figure_6.pdf", height = 4, width = 7)
par(mfrow = c(1,2), mar = c(5,4,4,1))

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate",
     main = "Individual hypotheses using\nnaive family-wise correction", asp = T)
for(i in 1:4) {
  lines(c(1,hyp_fpr_list_naive[[i]],0), c(1,hyp_tpr_list_naive[[i]],0), col = col_vec[i],
        lwd = lwd_vec[i])
}

legend("bottomright", c("0% flip", "10% flip", "25% flip", "75% flip"),
       bty="n", fill=col_vec)

plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate",
     main = "Individual hypotheses using\nnormalized statistic", asp = T)
for(i in 1:4) {
  lines(c(1,hyp_fpr_list[[i]],0), c(1,hyp_tpr_list[[i]],0), col = col_vec[i],
        lwd = lwd_vec[i])
}

legend("bottomright", c("0% flip", "10% flip", "25% flip", "75% flip"),
       bty="n", fill=col_vec)

graphics.off()
