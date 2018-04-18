source("../experiment/simulation_header.R")

trials <- 1000
vec <- rep(NA, trials)

library(huge)
set.seed(10)
d <- 40
L <- huge.generator(n = 50, d = d, graph = "hub", g = 4)
cov_mat <- L$sigma

for(i in 1:trials){
  set.seed(10*i)
  if(i %% floor(trials/10) == 0) cat('*')

  x <- MASS::mvrnorm(300, mu = rep(0,d), Sigma = cov_mat)
  y <- MASS::mvrnorm(300, mu = rep(0,d), Sigma = cov_mat)

  vec[i] <- longitudinalGM::cai_test(x, y, trials = 200, cores = 8)
}

save.image("../experiment/null_hypothesis.RData")

# plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
# lines(c(0,1),c(0,1), lwd = 2, col = "red")
