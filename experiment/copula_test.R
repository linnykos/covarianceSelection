rm(list=ls())
set.seed(10)
alpha <- 0.9; beta <- 1-alpha
d <- 10
mat <- matrix(beta, d, d)
mat[1:(d/2), 1:(d/2)] <- alpha
mat[(d/2+1):d, (d/2+1):d] <- alpha
diag(mat) <- 1
range(eigen(mat)$values)

x <- MASS::mvrnorm(5000, rep(0,d), mat)

# create the list of densities
load("../data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
idx <- sample(1:ncol(genexp), d)

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})

x2 <- nonparanormal_transformation(x, den_list, mean_vec = rep(0,d),
                                   sd_vec = rep(1,d))

# plot all the marginals
for(i in 1:d){
  hist(x2[,i], breaks = 100, col = "gray", main = i)
}

# plot the bivariate densities
combn_mat <- combn(d, 2)
for(i in 1:ncol(combn_mat)){
  plot(x2[,combn_mat[1,i]], x2[,combn_mat[2,i]], asp = T, pch = 16,
       col = rgb(0,0,0,0.5), main = paste0(combn_mat[,i], collapse = "-"))
}

cov(x2)
image(cov(x2))
