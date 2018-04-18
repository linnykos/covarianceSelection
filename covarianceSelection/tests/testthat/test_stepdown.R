context("Test stepdown")

## .compute_all_numerator_bootstrap is correct

test_that(".compute_all_numerator_bootstrap gives the same results as manual calculations", {
  set.seed(10)
  n <- 10; d <- 15; k <- 10
  dat_list <- lapply(1:k, function(x){matrix(rnorm(n*d),n,d)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})

  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)})
  idx <- which(lower.tri(diag(d), diag = T))
  remaining_idx <- 1:k

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, idx, remaining_idx)

  expect_true(all(sapply(res, length) == d*(d-1)/2+d))
  expect_true(all(sapply(res, is.numeric)))
  expect_true(!any(sapply(res, is.matrix)))
})

test_that(".compute_all_numerator_bootstrap gives the same results as manual calculations", {
  set.seed(10)
  n <- 15; d <- 10; k <- 5
  dat_list <- lapply(1:k, function(x){matrix(rnorm(n*d),n,d)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})

  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)})
  idx <- which(lower.tri(diag(d), diag = T))
  remaining_idx <- 1:k

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, idx, remaining_idx)

  res2 <- vector("list", k)
  dat_list <- lapply(dat_list, scale, center = T, scale = F)
  for(l in 1:k){
    mat <- matrix(0, d, d)
    for(i in 1:d){
      for(j in 1:i){
        mat[i,j] <- sum(noise_list[[l]] * (dat_list[[l]][,i] * dat_list[[l]][,j] - cov_list[[l]][i,j]))/n
      }
    }
    res2[[l]] <- mat[lower.tri(mat, diag = T)]
  }

  dif <- sum(sapply(1:k, function(x){sum(abs(res[[x]] - res2[[x]]))}))

  expect_true(dif < 1e-6)
})

test_that(".compute_all_numerator_bootstrap works with remaining_idx", {
  set.seed(10)
  n <- 25; d <- 10; k <- 5
  dat_list <- lapply(1:k, function(x){matrix(rnorm(n*d),n,d)})
  noise_list <- lapply(1:k, function(x){rnorm(n)})

  cov_list <- lapply(dat_list, function(x){(n-1)/n*stats::cov(x)})
  idx <- which(lower.tri(diag(d), diag = T))
  remaining_idx <- c(1,4,5)

  res <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, idx, remaining_idx)
  res2 <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, idx, 1:k)

  expect_true(length(res) == length(res2))
  for(i in c(1,4,5)){
    expect_true(sum(abs(res2[[i]] - res[[i]])) <= 1e-6)
  }
})

#######################

## .compute_all_test_stat is correct

test_that(".compute_all_test_stat works", {
  k <- 5
  num_list <- lapply(1:k, function(x){seq(x, x^2, length.out = 10)})
  denom_list <- lapply(1:k, function(x){x:(x+9)})

  res <- .compute_all_test_stat(num_list, denom_list)

  expect_true(length(res) == k*(k-1)/2)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
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
  dat_list <- lapply(1:k, function(x){matrix(rnorm(n*d),n,d)})
  cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)})
  denom_list <- .compute_all_denom(dat_list, cov_list)

  expect_true(all(sapply(denom_list, length) == d*(d-1)/2 + d))
  expect_true(all(sapply(denom_list, is.numeric)))
  expect_true(!any(sapply(denom_list, is.matrix)))
})

###############

## stepdown is correct

test_that("stepdown works", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 25)

  expect_true(is.numeric(res))
  expect_true(length(res) > 0)
  expect_true(length(res) <= 5*4/2)
})

test_that("stepdown works with no denominator", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 25, denominator = F)

  expect_true(is.numeric(res))
  expect_true(length(res) > 0)
  expect_true(length(res) <= 5*4/2)
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

  expect_true(sum(bool_vec[res]) >= length(res)-sum(bool_vec[res]))
  expect_true(sum(!bool_vec[-res]) >= length(bool_vec) - length(res) - sum(!bool_vec[-res]))
})

test_that("stepdown can reject with no denominator", {
  set.seed(10)
  dat_list <- lapply(1:6, function(x){
    if(x <= 3) matrix(rnorm(500), 50, 10) else matrix(rnorm(500,sd = 2), 50, 10)
  })
  res <- stepdown(dat_list, denominator = F, trials = 25)

  combn_mat <- combn(6, 2)
  bool_vec <- apply(combn_mat, 2, function(x){
    bool1 <- ifelse(x[1] <= 3, T, F)
    bool2 <- ifelse(x[2] <= 3, T, F)
    if(bool1 == bool2) T else F
  })

  expect_true(sum(bool_vec[res]) >= length(res)-sum(bool_vec[res]))
  expect_true(sum(!bool_vec[-res]) >= length(bool_vec) - length(res) - sum(!bool_vec[-res]))
})

test_that("stepdown rejects everything with alpha is 1", {
  set.seed(10)
  dat_list <- lapply(1:5, function(x){matrix(rnorm(100),10,10)})
  res <- stepdown(dat_list, trials = 25, alpha = 1)

  expect_true(length(res) == 0)
})

###################

## .shorten_combn is correct

test_that(".shorten_combn works", {
  mat <- combn(10,2)
  res <- .shorten_combn(mat, 5)

  expect_true(is.matrix(res))
  expect_true(nrow(res) == 2)
  expect_true(all(as.numeric(res) <= 5))

  bool_vec <- apply(mat, 2, function(x){if(all(x <= 5)) return(TRUE) else return(FALSE)})
  expect_true(all(sort(as.numeric(res)) == sort(as.numeric(mat[,which(bool_vec)]))))
})

test_that(".shorten_combn still returns a matrix even if only one column is selected", {
  mat <- combn(10,2)
  res <- .shorten_combn(mat, 2)

  expect_true(all(dim(res) == c(2,1)))
})
