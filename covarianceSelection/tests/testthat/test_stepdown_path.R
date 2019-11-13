context("Test stepdown path")

## stepdown_path is correct

test_that("stepdown_path works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown_path(dat_list, trials = 25, iterations = 20)

  expect_true(class(res) == "stepdown")
  expect_true(is.list(res))
  expect_true(all(names(res) == c("t_vec", "boot")))
  expect_true(length(res$t_vec) == 5*4/2)
  expect_true(is.numeric(res$t_vec))
  expect_true(!is.matrix(res$t_vec))

  expect_true(is.list(res$boot))
  expect_true(all(sapply(res$boot, is.numeric)))
  expect_true(all(sapply(res$boot, is.matrix)))
  expect_true(length(res$boot) == 20)

  dim_mat <- sapply(res$boot, dim)
  expect_true(all(dim_mat[1,] == 25))
  expect_true(all(dim_mat[2,] == 5*4/2))
})

test_that("stepdown_path works with probability", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown_path(dat_list, trials = 25, iterations = 20, prob = 0.5)
  
  expect_true(class(res) == "stepdown")
  expect_true(is.list(res))
  expect_true(all(names(res) == c("t_vec", "boot")))
  expect_true(length(res$t_vec) == 5*4/2)
  expect_true(is.numeric(res$t_vec))
  expect_true(!is.matrix(res$t_vec))
  
  expect_true(is.list(res$boot))
  expect_true(all(sapply(res$boot, is.numeric)))
  expect_true(all(sapply(res$boot, is.matrix)))
  expect_true(length(res$boot) == 20)
  
  dim_mat <- sapply(res$boot, dim)
  expect_true(all(dim_mat[1,] == 25))
  expect_true(all(dim_mat[2,] == 5*4/2))
})

test_that("stepdown_path works with squared = F", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown_path(dat_list, trials = 25, iterations = 20, prob = 0.5, squared = F)
  
  expect_true(class(res) == "stepdown")
  expect_true(is.list(res))
  expect_true(all(names(res) == c("t_vec", "boot")))
  expect_true(length(res$t_vec) == 5*4/2)
  expect_true(is.numeric(res$t_vec))
  expect_true(!is.matrix(res$t_vec))
  
  expect_true(is.list(res$boot))
  expect_true(all(sapply(res$boot, is.numeric)))
  expect_true(all(sapply(res$boot, is.matrix)))
  expect_true(length(res$boot) == 20)
  
  dim_mat <- sapply(res$boot, dim)
  expect_true(all(dim_mat[1,] == 25))
  expect_true(all(dim_mat[2,] == 5*4/2))
})

test_that("stepdown_path gives a different result when squared = T vs. squared = F", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res1 <- stepdown_path(dat_list, trials = 25, iterations = 20, prob = 0.5, squared = T)
  res2 <- stepdown_path(dat_list, trials = 25, iterations = 20, prob = 0.5, squared = F)
  
  expect_true(sum(abs(res1$t_vec - res2$t_vec)) > 0 )
})


#######################

## stepdown_choose is correct

test_that("stepdown_choose works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.05)

  expect_true(is.numeric(res$null_idx))
  expect_true(!is.matrix(res$null_idx))
  expect_true(all(res$null_idx > 0))
  expect_true(all(res$null_idx %% 1 == 0))
})

test_that("stepdown_choose works with return_pvalue", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.05, return_pvalue = T)
  
  expect_true(is.numeric(res$pval))
  expect_true(all(!is.na(res$pval)))
  expect_true(sum(res$pval) > 0)
  expect_true(length(res$pval) == 5*4/2)
})


test_that("stepdown_choose gets smaller as alpha gets bigger", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res1 <- stepdown_choose(obj, alpha = 0.05)
  res2 <- stepdown_choose(obj, alpha = 0.5)

  expect_true(all(res2$null_idx %in% res1$null_idx))
})

test_that("stepdown_choose can handle very large alpha", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.75)

  expect_true(length(res$null_idx) == 0)
})

test_that("stepdown_choose returns the same as stepdown", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  obj <- stepdown_path(dat_list, trials = 25, iterations = 20)
  res <- stepdown_choose(obj, alpha = 0.05)

  res2 <- stepdown(dat_list, trials = 25, alpha = 0.05)

  expect_true(all(sort(res$null_idx) == sort(res2$null_idx)))
})

