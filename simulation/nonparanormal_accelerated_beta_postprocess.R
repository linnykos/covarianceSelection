rm(list=ls())
load("../results/nonparanormal_accelerated_beta.RData")

median_mat <- matrix(NA, ncol = nrow(paramMat), nrow = 4)
for(i in 1:length(res)){
  tmp <- res[[i]][which(sapply(res[[i]], length) != 1)]
  
  our_vec <- sapply(tmp, function(x){x$error_our})
  base_vec <- sapply(tmp, function(x){x$error_base})
  all_vec <- sapply(tmp, function(x){x$error_all})
  oracle_vec <- sapply(tmp, function(x){x$error_oracle})
  
  median_mat[,i] <- c(mean(our_vec), mean(base_vec), mean(all_vec), mean(oracle_vec))
}

png("../figures/appendix_8b.png", height = 1400, width = 1400, res = 300, units ="px")
par(mar = c(5,4,4,1))

plot(NA, xlim = c(0,1), ylim = range(c(range(median_mat),85)), ylab = "Spectral error difference",
     xlab = expression(paste(beta, " level")), main = "Difference of spectral error\nversus Beta level")
lines(c(0,1), rep(median(median_mat[4,]),2), lty = 2, col = rgb(205,40,54, max = 255), lwd = 4)
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[1,], lwd = 4)
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[2,], col = rgb(106,164,248, max = 255), lwd = 4)
lines(seq(0,1,length.out = nrow(paramMat)), median_mat[3,], col = rgb(149,220,144, max = 255), lwd = 4)

legend("topleft", c("COBS", "Oracle", "Base", "All"),
       bty="n", fill= c("black", rgb(205,40,54, max = 255), rgb(106,164,248, max = 255), rgb(149,220,144, max = 255)))
graphics.off()