
rm(list=ls())
library(simulation)
library(covarianceSelection)

set.seed(10)
ncores <- 10
# doMC::registerDoMC(cores = ncores)
verbose <- F

trials <- 25
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 1000, seq(0, 1, length.out = 21), 0.9, 0.95))
colnames(paramMat) <- c("num_group1", "num_group2", "num_group3", "n", "d",
                        "percentage", "alpha", "gamma")

# collect all the marginal densities
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
set.seed(10)
idx <- sample(1:ncol(genexp), paramMat[1,"d"])

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})
rm(list = c("genexp"))

#############

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

generate_data <- function(covar_list, num_partition, n, den_list){
  k <- sum(num_partition)
  type_vec <- rep(1:3, times = num_partition)
  dat_list <- vector("list", k)
  d <- nrow(covar_list[[1]])
  
  dat_list <- lapply(1:k, function(i){ mvnfast::rmvn(n, rep(0, d), covar_list[[type_vec[i]]]) })
  
  # nonparanormal transform
  dat_list <- lapply(1:length(dat_list), function(x){
    covarianceSelection::nonparanormal_transformation(dat_list[[x]], den_list, 
                                                      mean_vec = rep(0, d),
                                                      sd_vec = sqrt(diag(covar_list[[type_vec[x]]])))
  })
  
  dat_list
}

################

rule <- function(vec){
  covar_list <- generate_covariance(d = vec["d"], percentage = vec["percentage"])
  
  if(verbose) print(paste0("Finish generating covariances: ", Sys.time()))
  
  dat <- generate_data(covar_list, num_partition = vec[1:3],  n = vec["n"], den_list)
  
  if(verbose) print(paste0("Finish generating data: ", Sys.time()))
  
  dat
}

criterion <- function(dat, vec, y){
  set.seed(10*y)
  dat2 <- rule(vec)
  idx_base <- c(1:3, vec["num_group1"]+1, vec["num_group1"] + vec["num_group2"] + 1)
  
  if(verbose) print(paste0("Running covariance tests: ", Sys.time()))
  
  set.seed(y)
  obj <- covarianceSelection::stepdown(dat, alpha = vec["alpha"], trials = 100, cores = ncores, verbose = F,
                                       squared = F)
  
  if(verbose) print(paste0("Starting paritions: ", Sys.time()))
  n <- sum(vec[1:3])
  g <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::add_edges(g, edges = combn_mat[,obj$null_idx])
  
  idx_our <-  covarianceSelection::clique_selection(g, threshold = vec["gamma"])[[1]]
  
  if(verbose) print(paste0("Forming datasets: ", Sys.time()))
  dat_our <- do.call(rbind, dat[idx_our])
  dat_base <- do.call(rbind, dat[idx_base])
  dat_oracle <- do.call(rbind, dat[c(1:vec["num_group1"])])
  dat2_oracle <- do.call(rbind, dat2[c(1:vec["num_group1"])])
  dat_all <- do.call(rbind, dat)
  
  if(verbose) print(paste0("Forming covariances: ", Sys.time()))
  cov_our <- stats::cov(dat_our)
  cov_base <- stats::cov(dat_base)
  cov_oracle <- stats::cov(dat_oracle)
  cov2_oracle <- stats::cov(dat2_oracle)
  cov_all <- stats::cov(dat_all)
  
  if(verbose) print(paste0("Computing errors: ", Sys.time()))
  error_our <- abs(mgcv::slanczos(cov_our - cov2_oracle, k = 1)$values[1])
  error_base <- abs(mgcv::slanczos(cov_base - cov2_oracle, k = 1)$values[1])
  error_oracle <- abs(mgcv::slanczos(cov_oracle - cov2_oracle, k = 1)$values[1])
  error_all <- abs(mgcv::slanczos(cov_all - cov2_oracle, k = 1)$values[1])
  
  list(partition_idx = obj$null_idx, idx_our = idx_our, error_our = error_our,
       error_base = error_base, error_oracle = error_oracle, error_all = error_all)
}

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/nonparanormal_accelerated_beta_tmp.RData",
                                        verbose = T)
save.image("../results/nonparanormal_accelerated_beta.RData")

