rm(list = ls())
load("../results/gaussian.RData")

# remove res that errored
for(i in 1:length(res)){
  rm_idx <- which(sapply(res[[i]], length) == 1)
  if(length(rm_idx) > 0) res[[i]] <- res[[i]][-rm_idx]
}

set.seed(10)
ncores <- 21
doMC::registerDoMC(cores = ncores)

func <- function(z){
  set.seed(10)
  
  n <- sum(paramMat[1,1:3])
  g <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::add_edges(g, edges = combn_mat[,z])
  tmp <- covarianceSelection::clique_selection(g, threshold = 0.95, verbose = F, time_limit = 60)
  
  if(length(res) > 1){
    len <- sapply(tmp, function(zz){length(intersect(zz, 1:paramMat[1,1]))})
    tmp <- tmp[[which.max(len)]]
  } else {
    tmp <- tmp[[1]]
  }
  tmp
}

## func(res[[1]][[1]]$indices_list[[5]])

partition_mat <- lapply(res, function(x){
  cat("\n")
  lapply(x, function(y){
    cat("*")
    foreach::"%dopar%"(foreach::foreach(i = 1:length(y$indices_list)), func(y$indices_list[[i]]))
  })
})

save.image("../results/gaussian_2.RData")

####################


rm(list=ls())
load("../results/gaussian_2.RData")

idx <- c(1:paramMat[1,1])
idx_all <- c(1:sum(paramMat[1,1:3]))

hyp_tpr_list <- lapply(partition_mat, function(x){ #loop over simulation settings
  sapply(x, function(y){ #loop over trials
    c(1, sapply(1:21, function(i){ #loop over alphas
      length(which(idx %in% y[[i]]))/length(idx)
    }), 0)
  })
})

hyp_fpr_list <- lapply(partition_mat, function(x){
  sapply(x, function(y){
    c(1, sapply(1:21, function(i){
      length(which(!y[[i]] %in% idx))/(length(idx_all) - length(idx))
    }), 0)
  })
})

.monotone <- function(x){
  n <- length(x)
  for(i in 1:(n-1)){
    if(x[i] > x[i+1]) x[i+1] = x[i]
  }
  x
}

# assign colors
colfunc <- colorRampPalette(c(rgb(205,40,54, max = 255), rgb(149,219,144, max = 255)))
col_vec <- colfunc(4)
lwd_vec <- c(5.5, 5, 4.5, 4)

png("../figures/appendix_6a.png", height = 1400, width = 1300, res = 300, units ="px")
par(mar = c(5,4,4,1))

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T, xlab = "False positive rate", ylab = "True positive rate",
     main = "Selected partitions via largest\nquasi-clique (RoC curve)")
for(k in 1:4){
  roc_our <- roc_region(hyp_tpr_list[[k]], hyp_fpr_list[[k]])
  lines(c(0, seq(0, 1, length.out = 21), 1), c(0, .monotone(roc_our$median_tpr), 1), col = col_vec[k], lwd = lwd_vec[k])
}
legend("bottomright", rev(c("0% level", "30% level", "60% level", "100% level")),
       bty="n", fill=rev(col_vec))
graphics.off()
