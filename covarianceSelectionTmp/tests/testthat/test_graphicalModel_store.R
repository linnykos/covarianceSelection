context("Test graphical model storing")

## .clean_glmnet is correct

test_that(".clean_glmnet works", {
  set.seed(10)
  x <- matrix(stats::rnorm(100*20),100,20)
  y <- stats::rnorm(100)
  fit <- glmnet::glmnet(x,y)

  res <- .clean_glmnet(fit)

  expect_true(!any(c("df", "dim", "dev.ratio", "nulldev", "npasses", "jerr", "call", "nobs")
              %in% names(res)))
  expect_true(inherits(res,"glmnet"))
})

test_that(".clean_glmnet output can still use coef", {
  set.seed(10)
  x <- matrix(stats::rnorm(100*20),100,20)
  y <- stats::rnorm(100)
  fit <- glmnet::glmnet(x,y)

  obj <- .clean_glmnet(fit)
  res <- as.numeric(stats::coef(obj, s = 0.05))

  expect_true(length(res) == 21)
  expect_true(is.numeric(res))
})

#####################

## .compute_reg_coefficients_save is correct

test_that(".compute_reg_coefficients_save works", {
  set.seed(10)
  dat <- matrix(stats::rnorm(100*5),100,5)

  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- .compute_reg_coefficients_save(dat, 1, filename_func)

  rm(list = c("dat", "tmp"))
  filename <- filename_func(1)
  load(filename)
  res <- as.numeric(stats::coef(res, s = 0.05))

  expect_true(length(res) == 5)
  expect_true(is.numeric(res))

  file.remove(filename)
})

##############

## graphicalModel_store is correct

test_that("graphicalModel_store works", {
  set.seed(10)
  dat <- matrix(stats::rnorm(100*5),100,5)

  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func)

  filename_vec <- sapply(1:5, filename_func)
  expect_true(all(filename_vec %in% list.files()))

  file.remove(filename_vec)
})

test_that("graphicalModel_store works with cross validation", {
  set.seed(10)
  dat <- matrix(stats::rnorm(100*5),100,5)

  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func, cv = T)

  filename_vec <- sapply(1:5, filename_func)
  expect_true(all(filename_vec %in% list.files()))

  file.remove(filename_vec)
})
