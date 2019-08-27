context("Test graphical model")

## graphicalModel is correct

test_that("graphicalModel works", {
  set.seed(10)
  dat <- matrix(rnorm(500), 50, 10)

  res <- graphicalModel(dat)

  expect_true(is.matrix(res))
  expect_true(ncol(res) == nrow(res))
  expect_true(sum(abs(res - t(res))) < 1e-5)
})

test_that("graphicalModel converges to identity", {
  set.seed(10)
  dat <- MASS::mvrnorm(1000, rep(0, 10), diag(10))

  res <- graphicalModel(dat)

  expect_true(sum(abs(res - diag(10)))/10^2 < 1e-2)
})

test_that("graphicalModel converges to something", {
  set.seed(10)
  L <- huge::huge.generator(n = 1000, d = 12, graph = "hub", g = 4, verbose = F)

  res <- graphicalModel(L$data)

  expect_true(sum(abs(res - L$omega))/12^2 < 0.1)
})

########################

## .compute_reg_coefficients_cv is correct

test_that(".compute_reg_coefficients_cv works", {
  set.seed(10)
  dat <- matrix(rnorm(500), 50, 10)

  res <- .compute_reg_coefficients_cv(dat)

  expect_true(is.list(res))
  expect_true(length(res) == ncol(dat))
  expect_true(all(sapply(res, length) == ncol(dat)))
})

#######################

## .compute_sigma_vec is correct

test_that(".compute_sigma_vec works", {
  set.seed(10)
  dat <- matrix(rnorm(500), 50, 10)

  coef_list <- .compute_reg_coefficients_cv(dat)
  coef_mat <- do.call(cbind, coef_list)
  res <- .compute_sigma_vec(dat, coef_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == ncol(dat))
})


