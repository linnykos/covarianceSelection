context("Test pairwise test numerically")

## pairwise_test is correct

test_that("pairwise_test gives uniform p values and puts the correct location", {
  len <- 5; n <- 50; d <- 3
  set.seed(10)

  dat.null <- vector("list", len)
  for(i in 1:len){
    dat.null[[i]] <- MASS::mvrnorm(n, rep(0,d), diag(d))
  }

  p.null <- pairwise_test(dat.null)

  dat.alt <- vector("list", len)
  for(i in 1:(len-1)){
    dat.alt[[i]] <- MASS::mvrnorm(n, rep(0,d), diag(d))
  }
  dat.alt[[len]] <- MASS::mvrnorm(n, rep(0,d), 5*diag(d))

  p.alt <- pairwise_test(dat.alt)

  for(i in 1:4){
    expect_true(p.alt[5,i] <= mean(p.alt[lower.tri(p.alt)]))
  }

  quant.vec <- c(0, 0.25, 0.5, 0.75, 1)
  unif.null <- sum(abs(quantile(p.null[lower.tri(p.null)], probs = quant.vec) - quant.vec))
  unif.alt <- sum(abs(quantile(p.alt[lower.tri(p.alt)], probs = quant.vec) - quant.vec))

  expect_true(unif.null < unif.alt)
})

test_that("pairwise_test works with two dat_lists", {
  len <- 3; n <- 50; d <- 3
  set.seed(10)

  dat_list <- vector("list", len)
  for(i in 1:len){
    dat_list[[i]] <- MASS::mvrnorm(n, rep(0,d), diag(d))
  }

  dat_list2 <- vector("list", len*2)
  for(i in 1:len){
    dat_list2[[i]] <- MASS::mvrnorm(n, rep(0,d), diag(d))
  }
  for(i in (len+1):(2*len)){
    dat_list2[[i]] <- MASS::mvrnorm(n, rep(0,d), 5*diag(d))
  }

  p_mat <- pairwise_test(dat_list, dat_list2)

  expect_true(all(dim(p_mat) == c(6,3)))

  quant_vec <- c(0, 0.25, 0.5, 0.75, 1)
  unif_null <- sum(abs(quantile(as.numeric(p_mat[1:3,]), probs = quant_vec) - quant_vec))
  unif_alt <- sum(abs(quantile(as.numeric(p_mat[4:6,]), probs = quant_vec) - quant_vec))
  expect_true(unif_null < unif_alt)
})
