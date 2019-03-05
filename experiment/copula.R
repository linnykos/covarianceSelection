rm(list=ls())
load("../data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)

den_list <- list(stats::density(genexp[,50]), stats::density(genexp[,100]))

set.seed(10)
x <- MASS::mvrnorm(5000, c(0,0), matrix(c(1,.5,.5,1),2,2))
plot(x[,1],x[,2], asp = T, pch = 16, col = rgb(0,0,0,0.2))

for(i in 1:ncol(x)){
  val <- stats::pnorm(x[,i], mean = 0, sd = 1)
  
  #make a quantile vector of the target density
  quantile_vec <- cumsum(den_list[[i]]$y)
  quantile_vec <- quantile_vec/max(quantile_vec)
  
  idx <- sapply(val, function(x){
    tmp <- x - quantile_vec
    tmp[tmp < 0] <- NA
    which.min(tmp)
  }) 
  
  difference <- (val - quantile_vec[idx])/(quantile_vec[idx+1]-quantile_vec[idx])
  
  x[,i] <- den_list[[i]]$x[idx] + difference*(den_list[[i]]$x[idx+1]-den_list[[i]]$x[idx])
}

plot(x[,1],x[,2], asp = T, pch = 16, col = rgb(0,0,0,0.2))
