context("Test Cai covariance test numerically")


test_that("The p-values are more uniform under the null", {
  trials <- 50
  d <- 3; n <- 50

  p_null <- numeric(trials)
  for(i in 1:trials){
    set.seed(10*i)
    x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
    y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
    p_null[i] <- cai_test(x,y, trials = 50)
  }

  p_alt <- numeric(trials)
  for(i in 1:trials){
    set.seed(10*i)
    x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
    y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = 5*diag(d))
    p_alt[i] <- cai_test(x,y, trials = 50)
  }

  unif_null <- sum(abs(quantile(p_null, probs = seq(0, 1, length.out = 11)) -
    seq(0, 1, length.out = 11)))
  unif_alt <- sum(abs(quantile(p_alt, probs = seq(0, 1, length.out = 11)) -
    seq(0, 1, length.out = 11)))
  expect_true(unif_null < unif_alt)
})
