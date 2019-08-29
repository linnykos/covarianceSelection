rm(list=ls())
library(microbenchmark)

p <- 3000; n <- 10

# testing .c_compute_sigma
set.seed(10)
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)

manual_cov <- function(dat){t(dat)%*%dat/nrow(dat)}

res <- microbenchmark(
  stats::cov(dat), manual_cov(dat), covarianceSelectionTmp:::c_compute_sigma(dat), times = 100
)
plot(res)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr      min       lq     mean
# stats::cov(dat) 172.1368 173.5394 175.4214
# manual_cov(dat) 176.6561 179.0858 180.8543
# covarianceSelectionTmp:::.c_compute_sigma(dat) 244.7507 247.3124 255.0377
# median       uq      max neval
# 174.2681 175.0944 252.8201   100
# 180.0477 181.9721 203.3202   100
# 249.3451 253.4832 324.8574   100

##############

# testing .c_compute_variance
rm(list=ls())
.compute_variance <- function(mat, cov_mat){
  n <- nrow(mat)
  
  mat <- scale(mat, center = TRUE, scale = FALSE)
  mat2 <- mat^2
  
  t(mat2)%*%mat2/n - cov_mat^2
}

set.seed(10)
p <- 3000; n <- 10
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- stats::cov(dat)

res <- microbenchmark(
  .compute_variance(dat, cov_mat), covarianceSelectionTmp:::c_compute_variance(dat, cov_mat), times = 100
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr      min       lq
# .compute_variance(dat, cov_mat) 250.8158 287.3850
# covarianceSelectionTmp:::.c_compute_variance(dat, cov_mat) 303.0551 333.0452
# mean   median       uq      max neval
# 301.1266 301.5825 315.0834 398.9428   100
# 401.5429 427.1848 443.3182 513.6221   100

######################

# testing .c_compute_bootSigma
rm(list=ls())
.compute_bootSigma <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

set.seed(10)
p <- 3000; n <- 10
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- stats::cov(dat)
noise_vec <- rnorm(nrow(dat))


res <- microbenchmark(
  .compute_bootSigma(dat, noise_vec, cov_mat), covarianceSelectionTmp:::c_compute_bootSigma(dat, noise_vec, cov_mat), times = 100
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr
# .compute_bootSigma(dat, noise_vec, cov_mat)
# covarianceSelectionTmp:::.c_compute_bootSigma(dat, noise_vec,      cov_mat)
# min       lq     mean   median       uq      max neval
# 206.7321 223.7780 253.7574 233.8668 288.6268 375.3291   100
# 294.9329 316.4361 355.8011 330.9611 368.9359 521.9608   100

#########################

# testing .c_compute_covStat
rm(list=ls())
.compute_covStat <- function(num_x, num_y, denom_x, denom_y, prob = 1){
  res <- (num_x - num_y)^2/(denom_x + denom_y)
  stats::quantile(abs(res), prob = prob)
}

set.seed(10)
p <- 6000; n <- 10
dat1 <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat1 <- scale(dat1, center = TRUE, scale = FALSE)
dat2 <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat2 <- scale(dat2, center = TRUE, scale = FALSE)
num1 <- covarianceSelection:::.c_compute_sigma(dat1); num2 <- covarianceSelection:::.c_compute_sigma(dat2)
dem1 <- covarianceSelection:::.c_compute_variance(dat1, num1); dem2 <- covarianceSelection:::.c_compute_variance(dat2, num2)

res <- microbenchmark(
  .compute_covStat(num1, num2, dem1, dem2), covarianceSelection:::.c_compute_covStat(num1, num2, dem1, dem2), times = 100
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: seconds
# expr      min
# .compute_covStat(num1, num2, dem1, dem2) 2.693438
# covarianceSelection:::.c_compute_covStat(num1, num2, dem1, dem2) 6.960124
# lq    mean   median       uq      max neval
# 2.791349 2.94196 2.834462 3.060822 4.071903   100
# 7.240974 7.49931 7.427881 7.743516 8.679723   100
