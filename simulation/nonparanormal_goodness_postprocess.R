rm(list=ls())
load("../results/nonparanormal_goodness.RData")

k <- 2
for(i in 1:3){
  vec <- sort(res[[k]][[1]][[i]])
  
  if(i == 2){
    tmp <- diff(vec)
    idx <- which(tmp > 1e-8)
    vec <- vec[idx]
  }
  
  plot(sort(vec), seq(0,1,length.out = length(vec)), asp = T, xlab = "Theoretical quantiles",
       pch = 16,
       ylab = "Observed quantiles", main = i)
  lines(c(0,1),c(0,1), col = "red", lty = 2)
}

