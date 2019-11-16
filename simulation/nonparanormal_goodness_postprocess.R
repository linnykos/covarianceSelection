rm(list=ls())
load("../results/nonparanormal_goodness.RData")


main_vec <- c("COBS", "Base", "All", "Oracle")
beta_vec <- c("0", "30", "60", "100")
for(k in 1:4){
  png(paste0("../figures/appendix_", k, ".png"), width = 3000, height = 800, res = 300, units = "px")
  par(mfrow = c(1,4), mar = c(4,4,4,1))
  for(i in 1:3){
    vec <- sort(res[[k]][[1]][[i]])
    if(i == 2){
      tmp <- diff(vec)
      idx <- which(tmp > 1e-8)
      vec <- vec[idx]
    }
    
    plot(sort(vec), seq(0,1,length.out = length(vec)), asp = T, xlab = "Theoretical quantiles",
         pch = 16,
         ylab = "Observed quantiles", main = paste0(main_vec[i], ": Beta ", beta_vec[k], "% level"))
    lines(c(0,1),c(0,1), col = "red", lty = 2)
  }
  
  vec <- sort(res[[1]][[1]][[4]])
  
  if(i == 2){
    tmp <- diff(vec)
    idx <- which(tmp > 1e-8)
    vec <- vec[idx]
  }
  
  plot(sort(vec), seq(0,1,length.out = length(vec)), asp = T, xlab = "Theoretical quantiles",
       pch = 16,
       ylab = "Observed quantiles", main = main_vec[4])
  lines(c(0,1),c(0,1), col = "red", lty = 2)
  
  graphics.off()
}




