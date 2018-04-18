rm(list = ls())
library(simulation)
library(covarianceSelection)
source("../simulation/simulation_helper.R")

generate_covariance <- function(d, percentage){
  covar_base <- .generate_block(d)
  covar_shuffle1 <- .shuffle(covar_base, percentage = percentage)
  covar_shuffle2 <- .shuffle(covar_base, percentage = percentage)

  list(covar_base = covar_base, covar_shuffle1 = covar_shuffle1,
       covar_shuffle2 = covar_shuffle2)
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

rule_closure <- function(){
  function(vec, ...){
    covar_list <- generate_covariance(d = vec["d"],
                                      percentage = vec["Shuffle.Percent"])

    generate_data(covar_list, num_partition = vec[1:3],  n = vec["n"])
  }
}

criterion_closure <- function(trials = 1000, cores = 10,
                              threshold = 0.95){
  function(dat, vec, y, ...){
    set.seed(y)
    res <- covarianceSelection::stepdown_path(dat, trials = trials, denominator = T,
                                         cores = cores, verbose = F)

    combn_mat <- utils::combn(sum(vec[1:3]), 2)
    alpha_vec <- seq(0, 1, length.out = 21)

    indices_list <- lapply(alpha_vec, function(alpha){
      covarianceSelection::stepdown_choose(res, alpha = alpha, verbose = F)
    })

    #our method
    partition_clique_list <- lapply(indices_list, function(x){
      g <- igraph::graph.empty(n = sum(vec[1:3]), directed = F)
      g <- igraph::add_edges(g, edges = combn_mat[,x])

      if(igraph::ecount(g) == 0) return(NA)
      res <- covarianceSelection::clique_selection(g, threshold = threshold,
                                              mode = "or", verbose = F,
                                              time_limit = 180)
      covarianceSelection::select_clique(res, 1:15, igraph::as_adj(g))
    })

    #spectral method
    partition_spectral_list <- lapply(indices_list, function(x){
      g <- igraph::graph.empty(n = sum(vec[1:3]), directed = F)
      g <- igraph::add_edges(g, edges = combn_mat[, x])
      covarianceSelection::spectral_selection(g, K = 3, target_idx = 1:15)
    })

    #naive pvalue method
    naive_pval_vec <- apply(combn_mat, 2, function(x){
      covarianceSelection::cai_test(dat[[x[1]]], dat[[x[2]]], trials = trials, cores = cores)
    })

    naive_indices_list <- lapply(alpha_vec, function(alpha){
      which(stats::p.adjust(naive_pval_vec, method = "bonferroni") >= alpha)
    })

    list(indices_list = indices_list, partition_clique_list = partition_clique_list,
         partition_spectral_list = partition_spectral_list,
         naive_pval_vec = naive_pval_vec, naive_indices_list = naive_indices_list)
  }
}


############

trials <- 20
paramMat <- as.matrix(expand.grid(15, 5, 5, 25, 50, c(0, 0.1, 0.25, 0.75)))
colnames(paramMat) <- c("Num.Group1", "Num.Group2", "Num.Group3", "n", "d",
                        "Shuffle.Percent")

rule <- rule_closure()
criterion <- criterion_closure()

res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                        as_list = T, filepath = "../results/RoC_tmp.RData")

save.image("../results/simulation_RoC.RData")


