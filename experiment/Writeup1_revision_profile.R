rm(list=ls())
set.seed(10)
library(microbenchmark)

# testing .c_compute_sigma
set.seed(10)
dat <- matrix(rnorm(15*200), nrow = 15, ncol = 200)
dat <- scale(dat, center = TRUE, scale = FALSE)

manual_cov <- function(dat){t(dat)%*%dat/nrow(dat)}

res <- microbenchmark(
  stats::cov(dat), manual_cov(dat), .c_compute_sigma(dat), times = 1000
)
res2 <- res; res2$time <- log(res2$time); plot(res2)


# > res
# Unit: microseconds
# expr     min       lq     mean   median       uq
# stats::cov(dat) 380.558 491.3825 696.2119 549.2375 652.5740
# manual_cov(dat) 428.790 538.3525 864.2238 586.5510 715.7915
# .c_compute_sigma(dat) 127.698 313.7290 524.1444 423.6790 493.4775
# max neval
# 24064.93  1000
# 27594.95  1000
# 25393.58  1000

##############

# testing .c_compute_variance

.compute_variance <- function(mat, cov_mat){
  n <- nrow(mat)
  
  mat <- scale(mat, center = TRUE, scale = FALSE)
  mat2 <- mat^2
  
  t(mat2)%*%mat2/n - cov_mat^2
}

set.seed(10)
dat <- matrix(rnorm(15*200), nrow = 15, ncol = 200)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- .c_compute_sigma(dat)

res <- microbenchmark(
  .compute_variance(dat, cov_mat), .c_compute_variance(dat, cov_mat), times = 1000
)

# > res
# Unit: microseconds
# expr     min       lq      mean
# .compute_variance(dat, cov_mat) 660.019 870.9345 2128.1334
# .c_compute_variance(dat, cov_mat) 147.137 331.7125  924.9418
# median       uq       max neval
# 1061.6565 1458.421 108866.08  1000
# 568.6095  712.437  40413.41  1000

######################

# testing .c_compute_bootSigma

.compute_bootSigma <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

set.seed(10)
dat <- matrix(rnorm(15*200), nrow = 15, ncol = 200)
dat <- scale(dat, center = TRUE, scale = FALSE)
cov_mat <- .c_compute_sigma(dat)
noise_vec <- rnorm(nrow(dat))


res <- microbenchmark(
  .compute_bootSigma(dat, noise_vec, cov_mat), .c_compute_bootSigma(dat, noise_vec, cov_mat), times = 1000
)

# > res
# Unit: microseconds
# expr     min       lq
# .compute_bootSigma(dat, noise_vec, cov_mat) 534.334 630.8165
# .c_compute_bootSigma(dat, noise_vec, cov_mat) 143.022 315.4490
# mean   median       uq      max neval
# 1106.2103 813.6655 972.5045 23437.23  1000
# 653.1823 491.5860 621.8845 21421.76  1000

#########################

# testing .c_compute_covStat

.compute_covStat <- function(num_x, num_y, denom_x, denom_y, prob = 1){
  res <- (num_x - num_y)^2/(denom_x + denom_y)
  stats::quantile(abs(res), prob = prob)
}

set.seed(10)
dat1 <- matrix(rnorm(15*200), nrow = 15, ncol = 200)
dat1 <- scale(dat1, center = TRUE, scale = FALSE)
dat2 <- matrix(rnorm(15*200), nrow = 15, ncol = 200)
dat2 <- scale(dat2, center = TRUE, scale = FALSE)
num1 <- .c_compute_sigma(dat1); num2 <- .c_compute_sigma(dat2)
dem1 <- .compute_variance(dat1, num1); dem2 <- .compute_variance(dat2, num2)

res <- microbenchmark(
  .compute_covStat(num1, num2, dem1, dem2), .c_compute_covStat(num1, num2, dem1, dem2), times = 1000
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: microseconds
# expr      min        lq
# .compute_covStat(num1, num2, dem1, dem2)  862.115  995.1685
# .c_compute_covStat(num1, num2, dem1, dem2) 2361.091 2649.6080
# mean   median       uq      max neval
# 2196.300 1255.619 1592.133 66908.06  1000
# 3566.146 3229.965 3752.394 32665.68  1000
