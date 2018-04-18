#plot the spectral difference in covariance
rm(list=ls())
load("../results/simulation_SnR.RData")

# select a particular alpha
idx <- which( c(0.1, 0.3, 0.7) == 0.7)
error_vec <- sapply(1:length(res), function(i){
  tmp <- sapply(1:trials, function(x){
    res[[i]][[x]]$error_vec[idx]
  })
  mean(tmp, na.rm = T)
})

# construct the curves
mat <- matrix(NA, ncol = 4, nrow = nrow(paramMat))
for(i in 1:length(res)){
  vec <- sapply(1:trials, function(x){
    c(res[[i]][[x]]$oracle_error, res[[i]][[x]]$base_error,
      res[[i]][[x]]$all_error)
  })
  mat[i,] <- c(error_vec[i], apply(vec, 1, mean))
}
mat <- mat[-nrow(mat),]
colnames(mat) <- c("Error", "Oracle.error", "Base.error", "All.error")

red <- rgb(205,40,54,maxColorValue=255)
blue <- rgb(106,144,202,maxColorValue=255)
green <- rgb(149,219,144,maxColorValue=255)
col_vec <- c("black", red, blue, green)
lty_vec <- c(1, 2, 1, 1)

pdf("../figures/figure_8.pdf", height = 4.25, width = 6.5)
ylim <- c(min(mat), max(mat))
plot(NA, xlim = range(paramMat[-nrow(paramMat),"Shuffle.Percent"])*100,
     ylim = ylim, main = "Difference of spectral error versus flip percentage",
     xlab = "Flip percentage", ylab = "Spectral error difference")
for(i in c(1:4)){
  lines(x = paramMat[-nrow(paramMat),"Shuffle.Percent"]*100, y = mat[,i],
        col = col_vec[i], lty = lty_vec[i], lwd = 5)
}

lines(c(-10, 10), y = rep(0, 2), lty = 2, lwd = 2)

legend("topleft", c("Stepdown", "Oracle", "Base", "All"),
       bty="n", fill=c("black", red, blue, green))
graphics.off()

# ################
#
# # find how much the best alpha varies
# seq_vec <- round(seq(0, 1, length.out = 21)*100)/100
# alpha_mat <- sapply(1:(length(res)-1), function(i){
#   sapply(1:trials, function(x){
#     seq_vec[which.min(unlist(res[[i]][[x]]$error_list))]
#   })
# })
#
# # compute the mean and sd for each level
# plot_info <- apply(alpha_mat, 2, function(x){
#   quantile(x, probs = c(0.25, 0.5, 0.75))
# })
#
# xvec <- paramMat[-nrow(paramMat),"Shuffle.Percent"]*100
# plot(NA, xlim = range(xvec), ylim = range(plot_info),
#      xlab = "Flip percentage", ylab = "Alpha level")
#
# lines(xvec, plot_info[2,], lwd = 5, col = red)
# for(i in 1:length(xvec)){
#   lines(rep(xvec[i], 2), plot_info[c(1,3),i], lwd = 3, col = red)
# }
