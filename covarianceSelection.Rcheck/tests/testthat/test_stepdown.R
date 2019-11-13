context("Test stepdown")

## .compute_all_numerator_bootstrap is correct

test_that(".compute_all_numerator_bootstrap gives the same results as manual calculations", {
  set.seed(10)
  n <- 10; d <- 15; k <- 10
  dat_list <- lapply(1:k, function(x){mat <- matrix(rnorm(n*d),n,d); scale(mat, center = T, scale = F)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))

  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)[diag_idx]})
  remaining_idx <- 1:k

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, diag_idx, remaining_idx)

  expect_true(all(sapply(res, length) == d*(d-1)/2+d))
  expect_true(all(sapply(res, is.numeric)))
})

test_that(".compute_all_numerator_bootstrap gives the same results as manual calculations", {
  set.seed(10)
  n <- 15; d <- 10; k <- 5
  dat_list <- lapply(1:k, function(x){mat <- matrix(rnorm(n*d),n,d); scale(mat, center = T, scale = F)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))

  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)})
  cov_list_trunc <- lapply(cov_list, function(x){x[diag_idx]})
  remaining_idx <- 1:k

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list_trunc, diag_idx, remaining_idx)

  res2 <- vector("list", k)
  for(l in 1:k){
    mat <- matrix(0, d, d)
    for(i in 1:d){
      for(j in 1:d){
        mat[i,j] <- sum(noise_list[[l]] * (dat_list[[l]][,i] * dat_list[[l]][,j] - cov_list[[l]][i,j]))/n
      }
    }
    res2[[l]] <- mat[diag_idx]
  }

  dif <- sum(sapply(1:k, function(x){sum(abs(res[[x]] - res2[[x]]))}))

  expect_true(dif < 1e-6)
})

test_that(".compute_all_numerator_bootstrap works with remaining_idx", {
  set.seed(10)
  n <- 25; d <- 10; k <- 5
  dat_list <- lapply(1:k, function(x){mat <- matrix(rnorm(n*d),n,d); scale(mat, center = T, scale = F)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  
  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)})
  cov_list_trunc <- lapply(cov_list, function(x){x[diag_idx]})
  idx <- which(lower.tri(diag(d), diag = T))
  remaining_idx <- c(1,4,5)

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list_trunc, diag_idx, remaining_idx)
  res2 <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list_trunc, diag_idx,  1:k)

  for(i in 1:3){
    expect_true(sum(abs(res2[[remaining_idx[i]]] - res[[i]])) <= 1e-6)
  }
})

#######################

## .compute_all_test_stat is correct

test_that(".compute_all_test_stat works", {
  set.seed(10)
  k <- 5
  num_list <- lapply(1:k, function(x){matrix(rnorm(25), 5, 5)})
  denom_list <- lapply(1:k, function(x){matrix(rnorm(25), 5, 5)})
  combn_mat <- combn(k, 2)

  res <- .compute_all_test_stat(num_list, denom_list, combn_mat)

  expect_true(length(res) == k*(k-1)/2)
})

test_that(".compute_all_test_stat changes with combn_mat", {
  k <- 5
  num_list <- lapply(1:k, function(x){seq(x, x^2, length.out = 10)})
  denom_list <- lapply(1:k, function(x){x:(x+9)})

  combn_mat <- combn(length(num_list), 2)
  combn_mat2 <- combn_mat[,c(1,2,5,8,10)]

  res <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat)
  res2 <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat2)

  expect_true(all(res[c(1,2,5,8,10)] == res2))
})

###############

## .compute_all_denom is correct

test_that(".compute_all_denom works", {
  set.seed(10)
  n <- 10; d <- 10; k <- 5
  dat_list <- lapply(1:k, function(x){mat <- matrix(rnorm(n*d),n,d); scale(mat, center = T, scale = F)})
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)[diag_idx]})
  
  res <- .compute_all_denom(dat_list, cov_list, diag_idx)

  expect_true(all(sapply(res, length) == d*(d-1)/2+d))
  expect_true(all(sapply(res, is.numeric)))
})

###############

## stepdown is correct

test_that("stepdown works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 25, return_pvalue = F)

  expect_true(is.numeric(res$null_idx))
  expect_true(length(res$null_idx) > 0)
  expect_true(length(res$null_idx) <= 5*4/2)
  expect_true(length(res$pval) == 1)
  expect_true(all(is.na(res$pval)))
})

test_that("stepdown can return the p values", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 100, return_pvalue = T)
  
  expect_true(is.numeric(res$null_idx))
  expect_true(length(res$null_idx) > 0)
  expect_true(length(res$null_idx) <= 5*4/2)
  expect_true(is.numeric(res$pval))
  expect_true(length(res$pval) == 5*4/2)
  expect_true(sum(res$pval) > 0)
  expect_true(all(!is.na(res$pval)))
})

test_that("stepdown can reject", {
  set.seed(10)
  dat_list <- lapply(1:6, function(x){
    if(x <= 3) matrix(rnorm(500), 50, 10) else matrix(rnorm(500,sd = 2), 50, 10)
  })
  res <- stepdown(dat_list, trials = 25)

  combn_mat <- combn(6, 2)
  bool_vec <- apply(combn_mat, 2, function(x){
    bool1 <- ifelse(x[1] <= 3, T, F)
    bool2 <- ifelse(x[2] <= 3, T, F)
    if(bool1 == bool2) T else F
  })

  expect_true(sum(bool_vec[res$null_idx]) >= length(res$null_idx)-sum(bool_vec[res$null_idx]))
  expect_true(sum(!bool_vec[-res$null_idx]) >= length(bool_vec) - length(res$null_idx) - sum(!bool_vec[-res$null_idx]))
})

test_that("stepdown rejects everything with alpha is 1", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 25, alpha = 1)

  expect_true(length(res$null_idx) == 0)
})
