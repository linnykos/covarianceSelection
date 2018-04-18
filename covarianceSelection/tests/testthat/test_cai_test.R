context("Test Cai covariance test")

## .compute_sigma is correct

test_that(".compute_sigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  res <- .compute_sigma(mat)

  mean.vec <- colSums(mat)/n
  res2.mat <- matrix(0, d, d)
  for(i in 1:n){
    res2.mat <- res2.mat + 1/n*(mat[i,] - mean.vec)%*%t(mat[i,] - mean.vec)
  }

  res2 <- res2.mat[lower.tri(res2.mat, diag = T)]

  expect_true(sum(abs(res - res2)) < 1e-4)
})

test_that(".compute_sigma gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  res <- .compute_sigma(mat)

  expect_true(length(res) == 5*4/2 + 5)
})

##############################

## .compute_variance is correct

test_that(".compute_variance is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  res <- .compute_variance(mat)

  sigma.mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2.mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:i){
      res2.mat[i,j] <- sum((mat[,i] * mat[,j] - sigma.mat[i,j])^2)/n
    }
  }

  res2 <- res2.mat[lower.tri(res2.mat, diag = T)]

  expect_true(sum(abs(res - res2)) < 1e-4)
})

test_that(".compute_variance is computed correctly when passing in cov and idx", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  cov_mat <- (n-1)/n*stats::cov(mat)
  idx <- which(lower.tri(diag(d), diag = T))
  res <- .compute_variance(mat, cov_mat, idx)
  res2 <- .compute_variance(mat)

  expect_true(sum(abs(res - res2)) < 1e-4)
})


test_that(".compute_variance is the same as an alternative way (in wrong order)", {
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
  res <- .compute_variance(mat)
  res2 <- func(mat)

  expect_true(sum(abs(sort(res) - sort(res2))) < 1e-4)
})

test_that(".compute_variance gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  res <- .compute_variance(mat)

  expect_true(length(res) == 5*4/2 + 5)
})
############################

## .compute_bootSigma is correct

test_that(".compute_bootSigma works", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  noise_vec <- rnorm(n)
  res <- .compute_bootSigma(mat, noise_vec)

  expect_true(is.numeric(res))
  expect_true(length(res) == d + d*(d-1)/2)
})

test_that(".compute_bootSigma is computed correctly", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  noise_vec <- rnorm(n)
  res <- .compute_bootSigma(mat, noise_vec)

  cov_mat <- (n-1)/n*stats::cov(mat)
  mat <- scale(mat, center = T, scale = F)
  res2.mat <- matrix(0, d, d)
  for(i in 1:d){
    for(j in 1:i){
      res2.mat[i,j] <- sum(noise_vec * (mat[,i] * mat[,j] - cov_mat[i,j]))/n
    }
  }

  res2 <- res2.mat[lower.tri(res2.mat, diag = T)]
  expect_true(sum(abs(res - res2)) < 1e-4)
})

test_that(".compute_variance gives the right length", {
  set.seed(10)
  d <- 5; n <- 10
  mat <- matrix(rnorm(n*d), n, d)
  res <- .compute_bootSigma(mat, rnorm(n))

  expect_true(length(res) == 5*4/2 + 5)
})

##############################

## .compute_covStat is correct

test_that(".compute_covStat returns 0 if nom.x == nom.y", {
  num_x <- 1:10; num_y <- 1:10
  denom_x <- 1:10; denom_y <- 1:10

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
