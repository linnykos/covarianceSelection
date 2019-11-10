rm(list=ls())
load("../results/nonparanormal_beta.RData")

median_mat <- matrix(NA, ncol = nrow(paramMat), nrow = 3)
for(i in 1:length(res)){
  our_vec <- sapply(res[[i]], function(x){x$error_our})
  base_vec <- sapply(res[[i]], function(x){x$error_base})
  all_vec <- sapply(res[[i]], function(x){x$error_all})
  
  median_mat[,i] <- c(median(our_vec), median(base_vec), median(all_vec))
}

plot(NA, xlim = c(0,1), ylim = range(c(0,range(median_mat))))
lines(c(0,1), rep(0,2), lty = 2, col = "red")
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[1,])
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[2,], col = "blue")
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[3,], col = "green")

