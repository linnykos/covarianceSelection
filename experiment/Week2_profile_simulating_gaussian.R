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

# > res
# Unit: seconds
# expr        min
# MASS::mvrnorm(n = 15, mu = rep(0, d), Sigma = covar_base) 121.918061
# mvnfast::rmvn(n = 15, mu = rep(0, d), sigma = covar_base)   7.627418
# lq       mean    median         uq        max neval
# 130.037732 133.326817 133.98518 135.746326 144.018445    10
# 7.699585   7.808839   7.76531   7.960578   8.019774    10