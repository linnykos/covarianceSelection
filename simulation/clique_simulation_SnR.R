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
    return(numeric(0))
  }
}

criterion_closure <- function(trials = 1000, cores = 10,
                              threshold = 0.95){
  function(dat, vec, y, ...){
    set.seed(y)
    covar_list <- generate_covariance(d = vec["d"],
                                      percentage = vec["Shuffle.Percent"])

    set.seed(y)
    dat_list <- generate_data(covar_list, num_partition = vec[1:3], n = vec["n"])

    alpha_vec <- c(0.1, 0.3, 0.7)

    partition_list <- vector("list", length = length(alpha_vec))
    error_vec <- numeric(length(alpha_vec))

    for(i in 1:length(alpha_vec)){
      index <- covarianceSelection::stepdown(dat_list, trials = trials, denominator = T,
                                        cores = cores, alpha = alpha_vec[i], verbose = F)

      edges <- utils::combn(length(dat_list), 2)

      g <- igraph::graph.empty(n = sum(vec[1:3]), directed = F)
      g <- igraph::add_edges(g, edges = edges[, index])

      if(igraph::ecount(g) == 0) next()
      res <- covarianceSelection::clique_selection(g, threshold = threshold,
                                              mode = "or", verbose = F,
                                              time_limit = 300)

      partition_list[[i]] <- covarianceSelection::select_clique(res, c(1:3,16,21),
                                                           igraph::as_adj(g))

      if(length(partition_list[[i]]) == 0) next()
      dat_list2 <- dat_list[partition_list[[i]]]
      dat_list2 <- lapply(dat_list2, scale, center = T, scale = F)
      dat_list2 <- do.call(rbind, dat_list2)
      cov_mat2 <- stats::cov(dat_list2)
      error_vec[i] <- .spectral_norm_mat(cov_mat2, covar_list$covar_base)
    }

    dat_list3 <- lapply(dat_list, scale, center = T, scale = F)
    dat_list3 <- do.call(rbind, dat_list3)
    cov_mat3 <- stats::cov(dat_list3)
    all_error <- .spectral_norm_mat(cov_mat3, covar_list$covar_base)

    dat_list4 <- dat_list[c(1:3,16,21)]
    dat_list4 <- lapply(dat_list4, scale, center = T, scale = F)
    dat_list4 <- do.call(rbind, dat_list4)
    cov_mat4 <- stats::cov(dat_list4)
    base_error <- .spectral_norm_mat(cov_mat4, covar_list$covar_base)

    dat_list5 <- dat_list[1:15]
    dat_list5 <- lapply(dat_list5, scale, center = T, scale = F)
    dat_list5 <- do.call(rbind, dat_list5)
    cov_mat5 <- stats::cov(dat_list5)
    oracle_error <- .spectral_norm_mat(cov_mat5, covar_list$covar_base)

    list(index = index,
         error_vec = error_vec, partition_list = partition_list,
         oracle_error = oracle_error, all_error = all_error,
         base_error = base_error)
  }
}

############

trials <- 10
paramMat <- as.matrix(expand.grid(15, 5, 5, 25, 50, seq(0, 1, length.out = 21)))
colnames(paramMat) <- c("Num.Group1", "Num.Group2", "Num.Group3", "n", "d",
                        "Shuffle.Percent")

rule <- rule_closure()
criterion <- criterion_closure()

res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                        as_list = T, filepath = "../results/SnR_tmp.RData")

save.image("../results/simulation_SnR.RData")

warnings()
