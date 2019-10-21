rm(list=ls())
library(simulation)
library(covarianceSelection)

paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 100, c(0, 0.3, 0.6, 1), 21))
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
  
  func <- function(i){
    mvnfast::rmvn(n, rep(0, d), covar_list[[typ_vec[i]]])
  }
  
  dat_list <- lapply(1:k, function(i){func(i)})
  
  # nonparanormal transform
  dat_list <- lapply(1:length(dat_list), function(x){
    covarianceSelection::nonparanormal_transformation(dat[[x]], den_list, 
                                                      mean_vec = rep(0, vec["d"]),
                                                      sd_vec = sqrt(diag(cov_list[[type_vec[i]]])))
  })
  
  dat_list
}

################

rule <- function(vec){
  covar_list <- generate_covariance(d = vec["d"], percentage = vec["percentage"])
  
  if(verbose) print(paste0("Finish generating covariances: ", Sys.time()))
  
  dat <- generate_data(covar_list, num_partition = vec[1:3],  n = vec["n"], den_list)
  
  if(verbose) print(paste0("Finish generating data: ", Sys.time()))
  
  dat_list
}

criterion <- function(dat, vec, y){
  set.seed(y)
  
  if(verbose) print(paste0("Starting to run the test: ", Sys.time()))
  
  obj <- covarianceSelection::stepdown_path(dat, trials = 200, cores = ncores, verbose = F,
                                            iterations = 10)
  tmp <- lapply(seq(0, 1, length.out = vec["alpha_levels"]), function(alpha){
    covarianceSelection::stepdown_choose(obj, alpha = alpha, return_pvalue = T)
  })
  
  # reformat
  indices_list <- vector("list", length =  vec["alpha_levels"]+1)
  for(i in 1:vec["alpha_levels"]){
    indices_list[[i]] <- tmp[[i]]$null_idx
  }
  names(indices_list) <- paste0("stepdown_", seq(0, 1, length.out = vec["alpha_levels"]))
  
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

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(2); criterion(rule(paramMat[10,]), paramMat[10,], 2)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = 20,
                                        cores = 15, as_list = T,
                                        filepath = "../results/nonparanormal_tmp.RData",
                                        verbose = T)
save.image("../results/nonparanormal.RData")