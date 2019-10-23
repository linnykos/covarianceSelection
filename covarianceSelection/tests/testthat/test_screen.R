context("Test screening")

## screen is correct

test_that("screen works", {
  set.seed(10)
  d <- 25
  dat <- huge::huge.generator(n = 50, d = d, graph = "hub", g = 4, verbose = F)$data
  pv <- c(stats::runif(10, max = 0.2), stats::runif(15))

  res <- screen(dat, pv, num_genes = 15)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(length(intersect(res$primary, res$secondary)) == 0)
})


test_that("screen outputs the correct number of variables", {
  set.seed(10)
  d <- 50
  dat <- huge::huge.generator(n = 50, d = d, graph = "hub", g = 4, verbose = F)$data
  pv <- c(stats::runif(10, max = 0.2), stats::runif(40))

  res <- screen(dat, pv, num_genes = 30)

  expect_true(length(res$primary) + length(res$secondary) == 30)
})

test_that("screen outputs the correct genes", {
  set.seed(10)
  cov_mat <- matrix(0, 20, 20)
  cov_mat[1:3, 1:3] <- 0.5
  cov_mat[4:20, 4:20] <- 0.5
  diag(cov_mat) <- 1
  pv <- c(0.01, 0.01, rep(0.5, 18))
  dat <- MASS::mvrnorm(100, mu = rep(0, 20), Sigma = cov_mat)

  res <- screen(dat, pv, num_genes = 3)

  expect_true(all(res$primary == c(1,2)))
  expect_true(res$secondary == 3)


})
