rm(list=ls())
# collect all the marginal densities
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
set.seed(10)
d <- 500
idx <- sample(1:ncol(genexp), d)

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})

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


set.seed(10)
covar_list <- generate_covariance(d = 500, percentage = 1)
dat <- generate_data(covar_list, num_partition = c(1,1,1),  n = 1000, den_list)
