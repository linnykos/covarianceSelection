rm(list=ls())
library(microbenchmark)
source("../simulation/simulation_helper.R")

d <- 3000
covar_base <- .generate_block(d, alpha = 0.9, beta = 0.1, spillover_percentage = 0,
                              normalize = F)

res <- microbenchmark::microbenchmark(
  MASS::mvrnorm(n = 15, mu = rep(0, d), Sigma = covar_base),
  mvnfast::rmvn(n = 15, mu = rep(0, d), sigma = covar_base), times = 10
)