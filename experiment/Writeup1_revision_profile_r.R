rm(list=ls())
library(microbenchmark)
.compute_bootSigma_1 <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

# .compute_bootSigma_2 <- function(mat, noise_vec, cov_mat){
#   n <- nrow(mat)
#   mat <- scale(mat, center = TRUE, scale = FALSE)
#   mat2 <- apply(mat, 2, function(x){x*noise_vec/n})
#   t(mat)%*%mat2 - (sum(noise_vec)/n)*cov_mat
# }
# 
# .compute_bootSigma_3 <- function(mat, noise_vec, cov_mat){
#   n <- nrow(mat)
#   mat <- scale(mat, center = TRUE, scale = FALSE)
#   mat2 <- apply(mat, 2, function(x){x*noise_vec/n})
#   crossprod(mat, mat2) - (sum(noise_vec)/n)*cov_mat
# }

.compute_bootSigma_2 <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  mat2 <- apply(mat, 2, function(x){x*noise_vec/n})
  Matrix::crossprod(mat, mat2) - (sum(noise_vec)/n)*cov_mat
}


set.seed(10)
p <- 1000; n <- 10
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- stats::cov(dat)
noise_vec <- rnorm(nrow(dat))

res1 <- .compute_bootSigma_1(dat, noise_vec, cov_mat)
res2 <- .compute_bootSigma_2(dat, noise_vec, cov_mat)
# res3 <- .compute_bootSigma_3(dat, noise_vec, cov_mat)

sum(abs(res1- res2))

res <- microbenchmark(
  .compute_bootSigma_1(dat, noise_vec, cov_mat), 
  .compute_bootSigma_2(dat, noise_vec, cov_mat), times = 100
)

# > res
# Unit: milliseconds
# expr      min
# .compute_bootSigma_1(dat, noise_vec, cov_mat) 8.647340
# .compute_bootSigma_2(dat, noise_vec, cov_mat) 9.682189
# lq     mean   median       uq       max neval
# 14.25275 19.35156 14.68832 19.73437  82.11067   100
# 15.54721 24.16936 16.60232 23.88548 235.80291   100

############################################

rm(list=ls())
.compute_covStat <- function(num_x, num_y, denom_x, denom_y, prob = 1){
  res <- (num_x - num_y)^2/(denom_x + denom_y)
  stats::quantile(abs(res), prob = prob)
}

.compute_covStat2 <- function(num_x, num_y, denom_x, denom_y, prob = 1){
  res <- (num_x - num_y)^2/(denom_x + denom_y)
  last_idx <- ceiling(prod(dim(num_x))*prob)
  sort(abs(res), partial = 1:last_idx)[last_idx]
}

set.seed(10)
p <- 500; n <- 10
dat1 <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat1 <- scale(dat1, center = TRUE, scale = FALSE)
dat2 <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat2 <- scale(dat2, center = TRUE, scale = FALSE)
num1 <- covarianceSelectionTmp:::c_compute_sigma(dat1); num2 <- covarianceSelectionTmp:::c_compute_sigma(dat2)
dem1 <- covarianceSelectionTmp:::c_compute_variance(dat1, num1); dem2 <- covarianceSelectionTmp:::c_compute_variance(dat2, num2)

.compute_covStat(num1, num2, dem1, dem2)
.compute_covStat2(num1, num2, dem1, dem2)

res <- microbenchmark::microbenchmark(
  .compute_covStat(num1, num2, dem1, dem2), 
  .compute_covStat2(num1, num2, dem1, dem2), times = 100
)
