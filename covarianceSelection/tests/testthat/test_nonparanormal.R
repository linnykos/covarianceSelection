context("Test nonparanormal")

test_that("nonparanormal_transformation works", {
  set.seed(10)
  tmp <- rexp(n = 10000)
  den_list <- list(stats::density(tmp))
  dat <- matrix(rnorm(15), ncol = 1)
  res <- nonparanormal_transformation(dat, den_list, mean_vec = 0, sd_vec = 1)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})
