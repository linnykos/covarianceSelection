rm(list=ls())
library(microbenchmark)

p <- 200; n <- 10

# testing .c_compute_sigma
set.seed(10)
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)

manual_cov <- function(dat){t(dat)%*%dat/nrow(dat)}

res <- microbenchmark(
  stats::cov(dat), manual_cov(dat), covarianceSelectionTmp:::.c_compute_sigma(dat), times = 1000
)
plot(res)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr      min        lq      mean
# stats::cov(dat) 625.7506  685.7971  715.4059
# manual_cov(dat) 679.3357  697.0153  727.2036
# covarianceSelection:::.c_compute_sigma(dat) 921.8879 1037.5564 1094.7013
# median        uq       max neval
# 696.8555  742.6572  819.5590   100
# 714.1300  748.7725  885.0109   100
# 1101.0343 1150.9371 1382.1910   100

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
p <- 6000; n <- 10
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- stats::cov(dat)

res <- microbenchmark(
  .compute_variance(dat, cov_mat), covarianceSelection:::.c_compute_variance(dat, cov_mat), times = 100
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr       min       lq
# .compute_variance(dat, cov_mat)  992.6049 1025.341
# covarianceSelection:::.c_compute_variance(dat, cov_mat) 1145.5202 1224.760
# mean   median       uq      max neval
# 1098.129 1062.467 1172.194 1370.378   100
# 1329.875 1300.299 1397.318 2030.893   100

######################

# testing .c_compute_bootSigma
rm(list=ls())
.compute_bootSigma <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

set.seed(10)
p <- 6000; n <- 10
dat <- matrix(rnorm(n*p), nrow = n, ncol = p)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- stats::cov(dat)
noise_vec <- rnorm(nrow(dat))


res <- microbenchmark(
  .compute_bootSigma(dat, noise_vec, cov_mat), covarianceSelection:::.c_compute_bootSigma(dat, noise_vec, cov_mat), times = 100
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: milliseconds
# expr       min
# .compute_bootSigma(dat, noise_vec, cov_mat)  850.4727
# covarianceSelection:::.c_compute_bootSigma(dat, noise_vec, cov_mat) 1219.6890
# lq      mean    median       uq      max neval
# 894.288  955.0559  954.4304 1003.158 1213.916   100
# 1393.235 1474.4309 1471.6795 1554.023 1834.537   100

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
