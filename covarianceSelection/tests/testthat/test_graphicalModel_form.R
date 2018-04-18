context("Test graphical model forming")

## .extract_beta is correct

test_that(".extract_beta works", {
  set.seed(10)
  d <- 12
  dat <- huge::huge.generator(n = 50, d = d, graph = "hub", g = 4, verbose = F)$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func)

  rm(list = c("dat", "tmp"))

  res <- .extract_beta(filename_func, 0.4, d)

  expect_true(is.list(res))
  expect_true(length(res) == d)
  expect_true(all(sapply(res, length) == d))
  expect_true(all(sapply(res, function(x){length(x == 0) >= 1})))
  expect_true(all(sapply(res, function(x){length(x != 0) >= 1})))

  filename_vec <- sapply(1:d, filename_func)
  file.remove(filename_vec)
})

test_that(".extract_beta works with num_primary", {
  set.seed(10)
  d <- 5
  num_primary <- 3
  gen <- huge::huge.generator(n = 500, d = d, graph = "hub", g = 4, verbose = F)
  dat <- gen$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func, num_primary = 3)

  res <- .extract_beta(filename_func, lambda = 0.01, d = d, num_primary = num_primary)

  expect_true(length(res) == 3)
  expect_true(all(sapply(res, length) == 5))
})

##################

## graphicalModel_form is correct

test_that("graphicalModel_form works", {
  set.seed(10)
  d <- 5
  gen <- huge::huge.generator(n = 500, d = d, graph = "hub", g = 4, verbose = F)
  dat <- gen$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func)

  res <- graphicalModel_form(dat, filename_func, 0.4)
  idx <- which(gen$omega == 0)
  expect_true(all(res[idx] == 0))

  filename_vec <- sapply(1:d, filename_func)
  file.remove(filename_vec)
})

test_that("graphicalModel_form works with cross validation", {
  set.seed(10)
  d <- 5
  gen <- huge::huge.generator(n = 500, d = d, graph = "hub", g = 4, verbose = F)
  dat <- gen$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func, cv = T)

  res <- graphicalModel_form(dat, filename_func)
  idx <- which(gen$omega == 0)
  expect_true(all(res[idx] == 0))

  filename_vec <- sapply(1:d, filename_func)
  file.remove(filename_vec)
})

test_that("graphicalModel_form works with num_primary", {
  set.seed(10)
  d <- 5
  num_primary <- 3
  gen <- huge::huge.generator(n = 500, d = d, graph = "hub", g = 4, verbose = F)
  dat <- gen$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func, num_primary = num_primary)

  res <- graphicalModel_form(dat, filename_func, 0.4, num_primary = num_primary)
  idx <- which(gen$omega == 0)
  expect_true(all(res[idx] == 0))

  expect_true(all(dim(res) == c(5,5)))
  expect_true(all(res == t(res)))

  filename_vec <- sapply(1:num_primary, filename_func)
  file.remove(filename_vec)
})

################

## .readjust_coef is correct

test_that(".readjust_coef works", {
  set.seed(10)
  d <- 5
  gen <- huge::huge.generator(n = 500, d = d, graph = "hub", g = 4, verbose = F)
  dat <- gen$data
  filename_func <- function(i){paste0("test_", i, ".RData")}
  tmp <- graphicalModel_store(dat, filename_func, num_primary = 3)

  coef_list <- .extract_beta(filename_func, lambda = 0.01, d = d, num_primary = 3)
  coef_mat <- do.call(cbind, coef_list)
  res <- .readjust_coef(coef_mat, d)

  expect_true(all(dim(res) == c(5,5)))
  expect_true(all(res[4:5,1:3] == t(res[1:3,4:5])))
})
