rm(list=ls())
source("../simulation/simulation_main.R")
library(longitudinalGM)

trials <- 1000
d <- 25; k <- 25; n <- 50
doMC::registerDoMC(cores = 14)

func1 <- function(i){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(i)
  covar_base <- .generate_block(d)
  dat_list <- lapply(1:k, function(y){
    MASS::mvrnorm(n, rep(0, d), covar_base)
  })

  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)})
  denom_list <- longitudinalGM:::.compute_all_denom(dat_list, cov_list)

  t_vec <- longitudinalGM:::.compute_all_test_stat(dat_list, denom_list, idx = diag_idx)

  max(t_vec)
}

sampling_dist <-as.numeric(unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                                                        func1(i))))

set.seed(1)
covar_base <- .generate_block(d)
dat_list <- lapply(1:k, function(y){
  MASS::mvrnorm(n, rep(0, d), covar_base)
})

diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)})
denom_list <- longitudinalGM:::.compute_all_denom(dat_list, cov_list)

len <- length(dat_list)
combn_mat <- utils::combn(len, 2)
idx_all <- rep(TRUE, ncol(combn_mat))

func2 <- function(i){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(i)
  noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
  remaining_pairs <- which(idx_all)
  combn_short <- combn_mat[, remaining_pairs, drop = F]
  if(any(is.na(combn_short))) return(Inf)

  remaining_idx <- unique(as.vector(combn_short))
  num_list <- longitudinalGM:::.compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, diag_idx,
                                               remaining_idx = remaining_idx)

  max(abs(longitudinalGM:::.compute_all_test_stat_bootstrap(num_list, denom_list, combn_mat = combn_short)))
}

bootstrap_dist <- as.numeric(unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                               func2(i))))

quantile(sampling_dist)
quantile(bootstrap_dist)

save.image("../experiment/bootstrap_dist.RData")
load("../experiment/bootstrap_dist.RData")

plot(sort(sampling_dist), sort(bootstrap_dist), asp = T)
lines(c(0,1), c(0,1), lwd = 2, col = "red")
