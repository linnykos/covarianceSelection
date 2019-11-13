context("Test goodness of fit")

## goodness_of_fit is correct

test_that("goodness_of_fit works", {
  set.seed(1)
  dat_list <- lapply(1:10, function(x){
    MASS::mvrnorm(50, mu = rep(0,4), Sigma = diag(4))
  })

  res <- goodness_of_fit(dat_list, 10, trials = 15)

  expect_true(length(res) == 10)
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})
