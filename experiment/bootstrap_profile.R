rm(list = ls())
library(simulation)
library(covarianceSelection)

paramMat <- cbind(5, 1, 1, 50, 100, c(0, 0.25, 0.5, 1))
colnames(paramMat) <- c("group1", "group2", "group3", "n", "d", "kappa")

# collect all the marginal densities
load("../data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
set.seed(10)
idx <- sample(1:ncol(genexp), paramMat[1,"d"])

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})

#############

generate_covariance1 <- function(vec){
  mat <- matrix(0.5, ncol = vec["d"], nrow = vec["d"])
  diag(mat) <- 1
  mat
}

generate_covariance2 <- function(vec){
  alpha <- 0.5+(0.95*vec["kappa"])*0.5
  beta <- 0.5-(0.95*vec["kappa"])*0.5
  mat <- matrix(beta, ncol = vec["d"], nrow = vec["d"])
  
  d <- vec["d"]
  d2 <- round(d/2)
  mat[1:d2, 1:d2] <- alpha
  mat[(d2+1):d, (d2+1):d] <- alpha
  diag(mat) <- 1
  
  mat
}

generate_covariance3 <- function(vec){
  generate_covariance1(vec)
}

find_cliques <- function(len, indices, threshold = 0.95){
  edges <- utils::combn(len, 2)
  
  g <- igraph::graph.empty(n = len, directed = F)
  g <- igraph::add_edges(g, edges = edges[,indices])
  
  if(igraph::ecount(g) == 0) return(NA)
  
  res <- covarianceSelection::clique_selection(g, threshold = threshold,
                                               mode = "or", verbose = F,
                                               time_limit = 300)
  
  covarianceSelection::select_clique(res, c(1:5,11,13), igraph::as_adj(g))
}

################

rule <- function(vec){
  covar1 <- generate_covariance1(vec)
  covar2 <- generate_covariance2(vec)
  covar3 <- generate_covariance3(vec)
  
  cov_list <- c(lapply(1:vec["group1"], function(x){covar1}), 
                lapply(1:vec["group2"], function(x){covar2}), 
                lapply(1:vec["group3"], function(x){covar3}))
  
  dat_list <- lapply(1:length(cov_list), function(x){
    MASS::mvrnorm(vec["n"], mu = rep(0, vec["d"]), Sigma = cov_list[[x]])
  })
  
  dat_list <- lapply(1:length(dat_list), function(x){
    covarianceSelection::nonparanormal_transformation(dat_list[[x]], den_list, 
                                                      mean_vec = rep(0, vec["d"]),
                                                      sd_vec = sqrt(diag(cov_list[[x]])))
  })
  
  for(i in (vec["group1"]+vec["group2"]+1):length(dat_list)){
    dat_list[[i]] <- (1+2*vec["kappa"])*dat_list[[i]]
  }
  
  dat_list
}

y <- 1
set.seed(y)
vec <- paramMat[2,]
dat <- rule(vec)

trials <- 1000
alpha_vec <- seq(0, 1, length.out = 21)
combn_mat <- utils::combn(sum(vec[1:3]), 2)

library(profvis)
p <- profvis({
  for(i in 1:10){
    set.seed(10*i)
    res <- covarianceSelection::stepdown(dat, trials = trials, 
                                         alpha = 0.05,
                                         denominator = T,
                                         cores = 1, verbose = F)
  }
})
save.image("../experiment/bootstrap_profile.RData")

htmlwidgets::saveWidget(p, paste0("../experiment/profile_large.html"))
# browseURL("../experiment/profile.html")