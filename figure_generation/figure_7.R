rm(list=ls())
load("../results/simulation_RoC.RData")

#########

# for partitions, compute TPR and FPR
# true positive rate
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

###########

graphsize <- rep(NA, 21)
clique_length <- rep(NA, 21)
spectral_length <- rep(NA, 21)
level <- 2
trial <- 2

for(i in 1:21){
  graphsize[i] <- length(res[[level]][[trial]]$indices_list[[i]])
  spectral_length[i] <- length(res[[level]][[trial]]$partition_spectral_list[[i]])
  clique_length[i] <- length(res[[level]][[trial]]$partition_clique_list[[i]])
}

############

pdf("../figures/figure_7.pdf", height = 4, width = 7)
par(mfrow = c(1,2), mar = c(5,4,4,1))

alpha_vec <- seq(0, 1, length.out = 21)
alpha_vec <- round(alpha_vec*100)/100
plot(NA, xlim = range(graphsize)-c(50,0), ylim = range(c(spectral_length, clique_length)),
     main = "Partitions selected based\non accepted hypothesis",
     ylab = "Number of selected partitions", xlab = "Number of accepted hypothesis")
points(x = graphsize, y = spectral_length, col = rgb(205,40,54, max = 255), cex = 1.25, lwd = 2)
points(x = graphsize, y = clique_length, col = "black", cex = 1.25, pch = 16)
legend("topleft", c("Largest partial clique", "Spectral clustering"),
       bty="n", fill= c("black", rgb(205,40,54, max = 255)))


plot(NA, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False positive rate", ylab = "True positive rate",
     main = "Selected partitions via largest\npartial clique (RoC curve)",
     asp = T)
for(j in 1:4){
  lines(c(1,par_fpr_list[[j]],0),
        c(1,par_tpr_list[[j]],0), col = col_vec[j], lwd = lwd_vec[j])
}


legend("bottomright", c("0% flip", "10% flip", "25% flip", "75% flip"),
       bty="n", fill=col_vec)

graphics.off()

