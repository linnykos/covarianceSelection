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
# stats::cov(dat) 160.9727 165.1366 171.5026
# manual_cov(dat) 167.2803 170.7199 177.1572
# covarianceSelectionTmp:::c_compute_sigma(dat) 247.6457 276.6433 293.4150
# median       uq      max neval
# 169.6042 175.5868 241.1468   100
# 175.5634 181.4187 204.6710   100
# 290.8825 302.0871 383.7218   100

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

# Unit: milliseconds
# expr      min       lq
# .compute_variance(dat, cov_mat) 250.4481 263.8595
# covarianceSelectionTmp:::c_compute_variance(dat, cov_mat) 246.7291 276.4481
# mean   median       uq      max neval
# 281.9083 271.1581 292.8610 418.8111   100
# 297.1633 295.7172 313.9755 451.8254   100

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

# res
# Unit: milliseconds
# expr
# .compute_bootSigma(dat, noise_vec, cov_mat)
# covarianceSelectionTmp:::c_compute_bootSigma(dat, noise_vec,      cov_mat)
# min       lq     mean   median       uq      max neval
# 208.1378 215.1320 231.7844 220.6916 246.7697 289.1676   100
# 247.4989 277.1277 291.1309 289.5158 309.6095 376.4795   100

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
num1 <- covarianceSelectionTmp:::c_compute_sigma(dat1); num2 <- covarianceSelectionTmp:::c_compute_sigma(dat2)
dem1 <- covarianceSelectionTmp:::c_compute_variance(dat1, num1); dem2 <- covarianceSelectionTmp:::c_compute_variance(dat2, num2)

res <- microbenchmark(
  .compute_covStat(num1, num2, dem1, dem2), covarianceSelectionTmp:::c_compute_covStat(num1, num2, dem1, dem2), times = 50
)
res2 <- res; res2$time <- log(res2$time); plot(res2)

# > res
# Unit: seconds
# expr
# .compute_covStat(num1, num2, dem1, dem2)
# covarianceSelectionTmp:::c_compute_covStat(num1, num2, dem1,      dem2)
# min       lq     mean   median       uq      max neval
# 2.772056 3.169486 3.291614 3.259967 3.447146 3.965732    50
# 6.121036 6.390369 6.593717 6.535867 6.854874 7.433616    50

#########################

rm(list=ls())
set.seed(10)
p <- 3000; n <- 10; r <- 10
dat_list <- lapply(1:r, function(x){matrix(rnorm(n*p), nrow = n, ncol = p)})
cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)})
noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
remaining_idx <- 1:length(dat_list)

.compute_bootSigma_tmp <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

.compute_all_numerator_bootstrap_tmp <- function(dat_list, noise_list, cov_list,
                                             remaining_idx){
  k <- length(dat_list)
  
  lis <- vector("list", k)
  lis[remaining_idx] <- lapply(remaining_idx, function(x){.compute_bootSigma_tmp(dat_list[[x]], noise_list[[x]],
                                                                             cov_list[[x]])})
  lis
}

# zz <- c_compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, remaining_idx)
res <- microbenchmark::microbenchmark(
  .compute_all_numerator_bootstrap_tmp(dat_list, noise_list, cov_list, remaining_idx), covarianceSelectionTmp:::c_compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, remaining_idx), times = 50
)
