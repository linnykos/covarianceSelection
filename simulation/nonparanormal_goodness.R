rm(list=ls())
library(simulation)
library(covarianceSelection)

set.seed(10)
ncores <- 20
doMC::registerDoMC(cores = ncores)
verbose <- F

trials <- 1
paramMat <- as.matrix(expand.grid(15, 5, 5, 15, 1000, c(0, 0.3, 0.6, 1), 0.1, 0.95))
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
  # dat_list <- lapply(1:length(dat_list), function(x){
  #   covarianceSelection::nonparanormal_transformation(dat_list[[x]], den_list, 
  #                                                     mean_vec = rep(0, d),
  #                                                     sd_vec = sqrt(diag(covar_list[[type_vec[x]]])))
  # })
  
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
  obj <- covarianceSelection::stepdown(dat, alpha = vec["alpha"], trials = 100, cores = ncores, verbose = F)
  n <- sum(vec[1:3])
  g <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::add_edges(g, edges = combn_mat[,obj$null_idx])
  idx_our <-  covarianceSelection::clique_selection(g, threshold = vec["gamma"], time_limit = 60)[[1]]
  
  goodness_our <- goodness_of_fit(dat[idx_our], permutations = 100, trials = 100, 
                                                       prob = 1, verbose = F)
  goodness_base <- covarianceSelection::goodness_of_fit(dat[c(1:3,16,21)], permutations = 250, trials = 100, 
                                                       prob = 1, verbose = F)
  goodness_all <- covarianceSelection::goodness_of_fit(dat, permutations = 250, trials = 100, 
                                                        prob = 1, verbose = F)
  
  if(vec["percentage"] == 0){
    goodness_oracle <- covarianceSelection::goodness_of_fit(dat[1:vec[1]], permutations = 250, trials = 100, 
                                                            prob = 1, verbose = F)
  } else {
    goodness_oracle <- NA
  }
  
  list(goodness_our = goodness_our, goodness_base = goodness_base,
       goodness_all = goodness_all, goodness_oracle = goodnes_oracle)
}
  
###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/nonparanormal_goodness_tmp.RData",
                                        verbose = T)
save.image("../results/nonparanormal_goodness.RData")
  
