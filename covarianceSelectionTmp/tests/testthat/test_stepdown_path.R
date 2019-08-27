context("Test stepdown path")

## stepdown_path is correct

test_that("stepdown_path works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown_path(dat_list, trials = 25, iterations = 20)

  expect_true(class(res) == "stepdown")
  expect_true(is.list(res))
  expect_true(all(names(res) == c("test", "boot")))
  expect_true(length(res$test) == 5*4/2)
  expect_true(is.numeric(res$test))
  expect_true(!is.matrix(res$test))

  expect_true(is.list(res$boot))
  expect_true(all(sapply(res$boot, is.numeric)))
  expect_true(all(sapply(res$boot, is.matrix)))
  expect_true(length(res$boot) == 20)

  dim_mat <- sapply(res$boot, dim)
  expect_true(all(dim_mat[1,] == 25))
  expect_true(all(dim_mat[2,] == 5*4/2))
})

test_that("stepdown_path works with no denominator", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown_path(dat_list, trials = 25, iterations = 20, denominator = F)

  expect_true(class(res) == "stepdown")
  expect_true(is.list(res))
  expect_true(all(names(res) == c("test", "boot")))
  expect_true(length(res$test) == 5*4/2)
  expect_true(is.numeric(res$test))
  expect_true(!is.matrix(res$test))

  expect_true(is.list(res$boot))
  expect_true(all(sapply(res$boot, is.numeric)))
  expect_true(all(sapply(res$boot, is.matrix)))
  expect_true(length(res$boot) == 20)

  dim_mat <- sapply(res$boot, dim)
  expect_true(all(dim_mat[1,] == 25))
  expect_true(all(dim_mat[2,] == 5*4/2))
})



#######################

## stepdown_choose is correct

test_that("stepdown_choose works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.05)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(all(res > 0))
  expect_true(all(res %% 1 == 0))
})

test_that("stepdown_choose gets smaller as alpha gets bigger", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res1 <- stepdown_choose(obj, alpha = 0.05)
  res2 <- stepdown_choose(obj, alpha = 0.5)

  expect_true(all(res2 %in% res1))
})

test_that("stepdown_choose can handle very large alpha", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.75)

  expect_true(length(res) == 0)
})

test_that("stepdown_choose returns the same as stepdown", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.05)

  res2 <- stepdown(dat_list, trials = 25, alpha = 0.05)

  expect_true(all(sort(res) == sort(res2)))
})

