rm(list = ls())
library(simulation)
library(covarianceSelection)
source("../simulation/simulation_helper.R")

trials <- 1
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 3000, 1))
colnames(paramMat) <- c("num_group1", "num_group2", "num_group3", "n", "d",
                        "percentage")

########3

generate_covariance <- function(d, percentage){
  covar_base <- .generate_block(d, alpha = 0.9, beta = 0.1, spillover_percentage = 0)
  covar_alt1 <- .generate_block(d, alpha = 0.9 - percentage*0.4, 
                                beta = 0.1 + percentage*0.4, 
                                spillover_percentage = 0)
  covar_alt2 <- .generate_block(d, alpha = 0.9, beta = 0.1, spillover_percentage = percentage*(1/6))
  
  list(covar_base = covar_base, covar_alt1 = covar_alt1,
       covar_alt2 = covar_alt2)
}

generate_data <- function(covar_list, num_partition, n){
  k <- sum(num_partition)
  type_vec <- rep(1:3, times = num_partition)
  dat_list <- vector("list", k)
  d <- nrow(covar_list[[1]])
  
  for(i in 1:k){
    if(type_vec[i] == 1) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_list[[1]])
    if(type_vec[i] == 2) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_list[[2]])
    if(type_vec[i] == 3) dat_list[[i]] <- MASS::mvrnorm(n, rep(0, d), covar_list[[3]])
  }
  
  dat_list
}

rule <- function(vec){
  covar_list <- generate_covariance(d = vec["d"], percentage = vec["percentage"])
  
  generate_data(covar_list, num_partition = vec[1:3],  n = vec["n"])

}

criterion <- function(dat, vec, y, ...){
  set.seed(y)
  res <- covarianceSelection::stepdown(dat, trials = 200, denominator = T, alpha = 0.5,
                                            cores = 15, verbose = T)
  
  list(res = res)
}

## set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

################

print(Sys.time())
res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                        as_list = T, filepath = "../experiment/test_tmp.RData")
print(Sys.time())




