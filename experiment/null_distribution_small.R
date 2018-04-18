source("../experiment/simulation_header.R")

trials <- 10
vec <- rep(NA, trials)

cov_mat <- diag(4)
cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

for(i in 1:trials){
  set.seed(10*i)
  if(i %% floor(trials/10) == 0) cat('*')

  x <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
  y <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)

  vec[i] <- longitudinalGM::cai_test(x,y)
}

save.image("../experiment/null_hypothesis_small.RData")

# plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
# lines(c(0,1),c(0,1), lwd = 2, col = "red")
