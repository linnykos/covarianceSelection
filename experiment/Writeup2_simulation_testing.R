rm(list = ls())
library(simulation)
library(covarianceSelection)
source("../simulation/simulation_helper.R")

set.seed(10)
ncores <- 15
doMC::registerDoMC(cores = ncores)

trials <- 10
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 1000, c(0, 0.5, 1), c(0, 0.5, 1)))
colnames(paramMat) <- c("num_group1", "num_group2", "num_group3", "n", "d",
                        "percentage", "alpha")

########3

generate_covariance <- function(d, percentage){
  covar_base <- .generate_block(d, alpha = 0.9, beta = 0.1, spillover_percentage = 0,
                                normalize = F)
  covar_alt1 <- .generate_block(d, alpha = 0.9 - percentage*0.4, 
                                beta = 0.1 + percentage*0.4, 
                                spillover_percentage = 0,
                                normalize = F)
  covar_alt2 <- .generate_block(d, alpha = 0.9, beta = 0.1, spillover_percentage = percentage*(1/6),
                                normalize = F)
  
  list(covar_base = covar_base, covar_alt1 = covar_alt1,
       covar_alt2 = covar_alt2)
}

generate_data <- function(covar_list, num_partition, n){
  k <- sum(num_partition)
  type_vec <- rep(1:3, times = num_partition)
  dat_list <- vector("list", k)
  d <- nrow(covar_list[[1]])
  
  func <- function(i){
    if(type_vec[i] == 1) return(mvnfast::rmvn(n, rep(0, d), covar_list[[1]]))
    if(type_vec[i] == 2) return(mvnfast::rmvn(n, rep(0, d), covar_list[[2]]))
    if(type_vec[i] == 3) return(mvnfast::rmvn(n, rep(0, d), covar_list[[3]]))
  }

  dat_list <- sapply(1:k, function(i){func(i)})
  
  # the nonparanormal transformation would happen here
  dat_list
}

rule <- function(vec){
  covar_list <- generate_covariance(d = vec["d"], percentage = vec["percentage"])
  
  print(paste0("Finish generating covariances: ", Sys.time()))
  
  dat <- generate_data(covar_list, num_partition = vec[1:3],  n = vec["n"])
  
  print(paste0("Finish generating data: ", Sys.time()))
  
  dat 
}

criterion <- function(dat, vec, y, ...){
  set.seed(y)
  
  print(paste0("Starting to run the test: ", Sys.time()))
  
  res <- covarianceSelection::stepdown(dat, trials = 200, denominator = T, alpha = vec["alpha"],
                                            cores = ncores, verbose = T)
  
  list(res = res)
}

## idx <- 8; y <- 3; set.seed(y); dat1 <- rule(paramMat[idx,]); set.seed(y); dat2 <- rule(paramMat[idx,])
## idx <- 8; y <- 3; set.seed(y); res <- criterion(rule(paramMat[idx,]), paramMat[idx,], y)

################

print(Sys.time())
res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                        as_list = T, filepath = "../experiment/test_tmp.RData")
print(Sys.time())

#################

# combn_null <- cbind(combn(paramMat[1,1],2),
#                     (combn(paramMat[1,2],2)+paramMat[1,1]),
#                     (combn(paramMat[1,3],2)+sum(paramMat[1,1:2])))
# num_partition <- sum(paramMat[1,1:3])
# idx_null <- combn_null[1,]+num_partition*combn_null[2,]
# combn_mat <- combn(num_partition,2)
# idx_all <- combn_mat[1,]+num_partition*combn_mat[2,]
# idx <- which(idx_all %in% idx_null)
# 
# z <- res[[1]][[1]]$res
# length(which(idx %in% z))/length(idx_null)
# length(which(1:paramMat[1,1] %in% z))/paramMat[1,1]


