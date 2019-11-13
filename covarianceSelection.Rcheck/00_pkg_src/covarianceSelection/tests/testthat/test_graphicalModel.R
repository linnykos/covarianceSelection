context("Test graphical model")

## graphicalModel is correct

test_that("graphicalModel works", {
  set.seed(10)
  dat <- matrix(rnorm(50^2), 50, 50)
  primary_idx <- 1:10
  lambda = 0.01

  res <- graphicalModel(dat, primary_idx, lambda)

  adj_mat <- as.matrix(res$adj_mat)
  expect_true(is.matrix(adj_mat))
  expect_true(ncol(adj_mat) == nrow(adj_mat))
  expect_true(sum(abs(adj_mat - t(adj_mat))) < 1e-5)
})

test_that("graphicalModel works with no specified lambda", {
  set.seed(10)
  dat <- matrix(rnorm(50^2), 50, 50)
  primary_idx <- 1:10
  
  res <- graphicalModel(dat, primary_idx, lambda = "lambda.1se")
  
  adj_mat <- as.matrix(res$adj_mat)
  expect_true(is.matrix(adj_mat))
  expect_true(ncol(adj_mat) == nrow(adj_mat))
  expect_true(sum(abs(adj_mat - t(adj_mat))) < 1e-5)
})

test_that("graphicalModel respects primary genes", {
  set.seed(10)
  dat <- matrix(rnorm(50^2), 50, 50)
  primary_idx <- 1:10
  lambda = 0.01
  
  res <- graphicalModel(dat, primary_idx, lambda)
  
  adj_mat <- as.matrix(res$adj_mat)
  expect_true(all(sum(abs(adj_mat[11:50, 11:50])) <= 1e-6))
  expect_true(sum(abs(adj_mat)) > 1e-6)
})
