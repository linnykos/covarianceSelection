rm(list=ls())
source("../simulation/simulation_main.R")
library(longitudinalGM)

generate_rule_closure <- function(n = 1000, d = 50){
  function(paramVec){
    k <- sum(paramVec[1:3])
    type_vec <- rep(1:3, times = paramVec[1:3])

    dat_list <- vector("list", k)
    covar_base <- .generate_block(d)
    covar_shuffle1 <- .shuffle(covar_base, percentage = paramVec[4])
    covar_shuffle2 <- .shuffle(covar_base, percentage = paramVec[4])
    for(i in 1:k){
      if(type_vec[i] == 1) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_base)
      if(type_vec[i] == 2) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_shuffle1)
      if(type_vec[i] == 3) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_shuffle2)
    }

    dat_list
  }
}

rule <- generate_rule_closure(n = 1000, d = 50)
paramVec <- c(15,5,5,0)
set.seed(0)
dat_list <- rule(paramVec)
res <- longitudinalGM::stepdown(dat_list, alpha = 1, trials = 10000, cores = 14, verbose = T)

save.image("../experiment/low_alpha.RData")
