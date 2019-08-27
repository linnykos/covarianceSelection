context("Test Cai covariance test")

## .c_compute_sigma is correct

test_that(".c_compute_sigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  res <- .c_compute_sigma(mat)

  mean.vec <- colSums(mat)/n
  res2.mat <- matrix(0, d, d)
  for(i in 1:n){
    res2.mat <- res2.mat + 1/n*(mat[i,] - mean.vec)%*%t(mat[i,] - mean.vec)
  }
  
  expect_true(sum(abs(res - res2.mat)) < 1e-4)
})

test_that(".c_compute_sigma gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  res <- .c_compute_sigma(mat)

  expect_true(length(res) == 5*5)
})

##############################

## .c_compute_variance is correct

test_that(".c_compute_variance is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  res <- .c_compute_variance(mat, cov_mat)

  sigma.mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2.mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:d){
      res2.mat[i,j] <- sum((mat[,i] * mat[,j] - sigma.mat[i,j])^2)/n
    }
  }

  expect_true(sum(abs(res - res2.mat)) < 1e-4)
})

test_that(".c_compute_variance is the same as an alternative way (in wrong order)", {
  func <- function(mat){
    d <- ncol(mat); n <- nrow(mat)
    mat <- scale(mat, center = TRUE, scale = FALSE)

    comb.mat <- utils::combn(d, 2)
    comb.mat <- cbind(comb.mat, matrix(rep(1:d, each = 2), ncol = d))

    s.vec <- apply(comb.mat, 2, function(x){
      zij <- mat[,x[1]] * mat[,x[2]]
      (n-1)/n * stats::var(zij)
    })

    s.vec
  }

  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  res <- .c_compute_variance(mat, cov_mat)
  res <- res[lower.tri(res, diag = T)]
  res2 <- func(mat)

  expect_true(sum(abs(sort(res) - sort(res2))) < 1e-4)
})

test_that(".c_compute_variance gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  res <- .c_compute_variance(mat, cov_mat)

  expect_true(length(res) == 5*5)
})

############################

## .c_compute_bootSigma is correct

test_that(".c_compute_bootSigma works", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  noise_vec <- rnorm(n)
  res <- .c_compute_bootSigma(mat, noise_vec, cov_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == d^2)
})

test_that(".c_compute_bootSigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  noise_vec <- rnorm(n)
  res <- .compute_bootSigma(mat, noise_vec, cov_mat)

  cov_mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2.mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:i){
      res2.mat[i,j] <- sum(noise_vec * (mat[,i] * mat[,j] - cov_mat[i,j]))/n
    }
  }
  res2 <- res2.mat[lower.tri(res2.mat, diag = T)]

  expect_true(sum(abs(sort(res) - sort(res2))) < 1e-4)
})

test_that(".c_compute_bootSigma gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- apply(mat, 2, function(x){x-mean(x)})
  cov_mat <- .c_compute_sigma(mat)
  res <- .compute_bootSigma(mat, rnorm(n), cov_mat)

  expect_true(length(res) == 5*4/2 + 5)
})

##############################

## .c_compute_covStat is correct

test_that("c_compute_covStat returns 0 if num_x == num_y", {
  num_x <- matrix(1:25,5,5); num_y <- matrix(1:25,5,5)
  denom_x <-  matrix(1:25,5,5); denom_y <- matrix(1:25,5,5)

  expect_true(.c_compute_covStat(num_x, num_y, denom_x, denom_y) == 0)
})

########################

## cai_test is correct

test_that("cai_test works", {
  d <- 5
  n <- 20
  set.seed(10)
  x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  res <- cai_test(x,y, trials = 50)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
  expect_true(res <= 1)
})

test_that("cai_test works with cores", {
  d <- 5
  n <- 20
  set.seed(10)
  x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  res <- cai_test(x,y, trials = 50, cores = 2)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
  expect_true(res <= 1)
})
