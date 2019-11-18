rm(list=ls())
library(simulation)
library(covarianceSelection)

set.seed(10)
ncores <- 10
doMC::registerDoMC(cores = ncores)
verbose <- F

trials <- 25
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 1000, c(0, 0.3, 0.6, 1), 21))
colnames(paramMat) <- c("num_group1", "num_group2", "num_group3", "n", "d",
                        "percentage", "alpha_levels")

# collect all the marginal densities
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
set.seed(10)
idx <- sample(1:ncol(genexp), paramMat[1,"d"])

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})
rm(list = "genexp")

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
  set.seed(y)
  alpha_vec <- seq(0, 1, length.out = vec["alpha_levels"])
  
  if(verbose) print(paste0("Starting to run the test: ", Sys.time()))
  
  obj <- covarianceSelection::stepdown_path(dat, trials = 200, cores = ncores, verbose = F,
                                            iterations = 10)
  tmp <- lapply(alpha_vec, function(alpha){
    covarianceSelection::stepdown_choose(obj, alpha = alpha, return_pvalue = T)
  })
  
  # reformat
  indices_list <- vector("list", length = length(alpha_vec))
  for(i in 1:length(alpha_vec)){
    indices_list[[i]] <- tmp[[i]]$null_idx
  }
  names(indices_list) <- paste0("stepdown_", seq(0, 1, length.out = length(alpha_vec)))
  
  naive_pval_vec <- tmp[[1]]$pval
  
  bonferroni_indices_list <- lapply(alpha_vec, function(alpha){
    which(stats::p.adjust(naive_pval_vec, method = "bonferroni") >= alpha)
  })
  
  bh_indices_list <- lapply(alpha_vec, function(alpha){
    which(stats::p.adjust(naive_pval_vec, method = "BH") >= alpha)
  })
  
  list(indices_list = indices_list, naive_pval_vec = naive_pval_vec,
       bonferroni_indices_list = bonferroni_indices_list,
       bh_indices_list = bh_indices_list)
}

# idx <- 1; y <- 1; set.seed(y); criterion(rule(paramMat[idx,]), paramMat[idx,], y)
# set.seed(2); criterion(rule(paramMat[10,]), paramMat[10,], 2)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/nonparanormal_tmp.RData",
                                        verbose = T)
save.image("../results/nonparanormal.RData")