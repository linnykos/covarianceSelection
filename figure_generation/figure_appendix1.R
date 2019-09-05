rm(list=ls())
load("../results/simulation_RoC_nodenom.RData")

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

# true postive rate
idx <- which(idx_all %in% idx_null)
hyp_tpr_list <- lapply(res, function(x){ #on each flip-rate
  mat <- sapply(x, function(y){ #trials
    sapply(y$indices_list, function(z){ #over alpha
      length(which(idx %in% z))/length(idx_null)
    })
  })

  apply(mat, 1, mean)
})

# false positive rate
hyp_fpr_list <- lapply(res, function(x){
  mat <- sapply(x, function(y){
    sapply(y$indices_list, function(z){
      length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
    })
  })

  apply(mat, 1, mean)
})

########

par_tpr_list <- lapply(res, function(x){ #on each flip-rate
  mat <- sapply(x, function(y){ #trials
    sapply(y$partition_clique_list, function(z){ #over alpha
      length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]
    })
  })

  apply(mat, 1, mean)
})

# false positive rate
par_fpr_list <- lapply(res, function(x){
  mat <- sapply(x, function(y){
    sapply(y$indices_list, function(z){
      length(which((paramMat[1,1]+1):(sum(paramMat[1,1:3])) %in% z))/(sum(paramMat[1,1:3])-paramMat[1,1])
    })
  })

  apply(mat, 1, mean)
})

###########

# create color ramp
colfunc <- colorRampPalette(c(rgb(205,40,54, max = 255), rgb(149,219,144, max = 255)))
col_vec <- colfunc(4)
lwd_vec <- c(5, 4.5, 4, 3.5)

pdf("../figures/figure_appendix1.pdf", height = 5, width = 5)
#par(mfrow = c(1,2), mar = c(5,4,4,1))
plot(NA, xlim = c(0,1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate",
     main = "Individual hypotheses using\nnon-normalized statistic", asp = T,
     cex.axis = 1.25, cex.lab = 1.25)
for(i in 1:4) {
  lines(c(1,hyp_fpr_list[[i]],0), c(1,hyp_tpr_list[[i]],0), col = col_vec[i],
        lwd = lwd_vec[i])
}

legend("bottomright", c("0% flip", "10% flip", "25% flip", "75% flip"),
       bty="n", fill=col_vec)

graphics.off()
