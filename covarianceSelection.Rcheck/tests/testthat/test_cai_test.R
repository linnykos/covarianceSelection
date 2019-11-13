context("Test Cai covariance test")
doMC::registerDoMC(cores = 1)

## .compute_sigma is correct

test_that(".compute_sigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  res <- .compute_sigma(mat, diag_idx)

  mean_vec <- colSums(mat)/n
  res2_mat <- matrix(0, d, d)
  for(i in 1:n){
    res2_mat <- res2_mat + 1/n*(mat[i,] - mean_vec)%*%t(mat[i,] - mean_vec)
  }
  
  expect_true(sum(abs(res - res2_mat[diag_idx])) < 1e-4)
})

test_that(".compute_sigma gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  res <- .compute_sigma(mat, diag_idx)

  expect_true(length(res) == length(diag_idx))
})

##############################

## .compute_variance is correct

test_that(".compute_variance is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  cov_mat <- .compute_sigma(mat, diag_idx)
  res <- .compute_variance(mat, cov_mat, diag_idx)

  sigma_mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2_mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:d){
      res2_mat[i,j] <- sum((mat[,i] * mat[,j] - sigma_mat[i,j])^2)/n
    }
  }

  expect_true(sum(abs(res - res2_mat[diag_idx])) < 1e-4)
})

test_that(".compute_variance is the same as an alternative way (in wrong order)", {
  func <- function(mat){
    d <- ncol(mat); n <- nrow(mat)
    mat <- scale(mat, center = TRUE, scale = FALSE)

    comb_mat <- utils::combn(d, 2)
    comb_mat <- cbind(comb_mat, matrix(rep(1:d, each = 2), ncol = d))

    s_vec <- apply(comb_mat, 2, function(x){
      zij <- mat[,x[1]] * mat[,x[2]]
      (n-1)/n * stats::var(zij)
    })

    s_vec
  }

  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  
  cov_mat <- .compute_sigma(mat, diag_idx)
  res <- .compute_variance(mat, cov_mat, diag_idx)
  res2 <- func(mat)

  expect_true(sum(abs(sort(res) - sort(res2))) < 1e-4)
})

test_that(".compute_variance gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  
  cov_mat <- .compute_sigma(mat, diag_idx)
  res <- .compute_variance(mat, cov_mat, diag_idx)

  expect_true(length(res) == length(diag_idx))
})

############################

## .compute_bootSigma is correct

test_that(".compute_bootSigma works", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  
  cov_mat <- .compute_sigma(mat, diag_idx)
  noise_vec <- rnorm(n)
  res <- .compute_bootSigma(mat, noise_vec, cov_mat, diag_idx)

  expect_true(is.numeric(res))
  expect_true(length(res) == length(diag_idx))
})

test_that(".compute_bootSigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  
  cov_mat <- .compute_sigma(mat, diag_idx)
  noise_vec <- rnorm(n)
  res <- .compute_bootSigma(mat, noise_vec, cov_mat, diag_idx)

  cov_mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2_mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:i){
      res2_mat[i,j] <- sum(noise_vec * (mat[,i] * mat[,j] - cov_mat[i,j]))/n
    }
  }
  res2 <- res2_mat[lower.tri(res2_mat, diag = T)]

  expect_true(sum(abs(sort(res) - sort(res2))) < 1e-4)
})

test_that(".compute_bootSigma gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  mat <- scale(mat, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(mat)), diag = T))
  
  cov_mat <- .compute_sigma(mat, diag_idx)
  res <- .compute_bootSigma(mat, rnorm(n), cov_mat, diag_idx)

  expect_true(length(res) == length(diag_idx))
})

##############################

## .compute_covStat is correct

test_that(".compute_covStat returns 0 if num_x == num_y", {
  num_x <- matrix(1:25,5,5); num_y <- matrix(1:25,5,5)
  denom_x <-  matrix(1:25,5,5); denom_y <- matrix(1:25,5,5)

  expect_true(.compute_covStat(num_x, num_y, denom_x, denom_y) == 0)
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