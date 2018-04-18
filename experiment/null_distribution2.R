rm(list=ls())
library(longitudinalGM)

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

d <- 50
n <- 25
trials <- 2
partition_vec <- c(15,5,5)

covar_base <- .generate_block(d)
pval_vec_null <- vector("list", length = trials)
pval_vec_alt <- vector("list", length = trials)

# for(i in 1:trials){
#   print(i)
#   set.seed(i)
#   covar_list <- generate_covariance(d = d, percentage = 0)
#   dat_list <- generate_data(covar_list, num_partition = partition_vec, n = n)
#
#   combn_mat <- utils::combn(sum(partition_vec), 2)
#
#   pval_vec_null[[i]] <- apply(combn_mat, 2, function(x){
#     longitudinalGM::cai_test(dat_list[[x[1]]], dat_list[[x[2]]], trials = 1000, cores = 10,
#                              prob = 0.9975)
#   })
# }

for(i in 1:trials){
  print(i)
  set.seed(i)
  covar_list <- generate_covariance(d = d, percentage = 0.7)
  dat_list <- generate_data(covar_list, num_partition = partition_vec, n = n)

  combn_mat <- utils::combn(sum(partition_vec), 2)

  pval_vec_alt[[i]] <- apply(combn_mat, 2, function(x){
    longitudinalGM::cai_test(dat_list[[x[1]]], dat_list[[x[2]]], trials = 1000, cores = 10,
                             prob = 1)
  })
}

save.image("../experiment/null_distribution2.RData")

load("../experiment/null_distribution2.RData")
# par(mfrow = c(1,2))
# vec <- unlist(pval_vec_null)
# plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
# lines(c(0,1), c(0,1), col = "red")

vec <- unlist(pval_vec_alt)
plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
lines(c(0,1), c(0,1), col = "red")

##############

#now plot the ROC curve
num_partition <- sum(partition_vec)

# for indices, compute TPR and FPR
combn_null <- cbind(combn(partition_vec[1],2),
                    (combn(partition_vec[2],2)+partition_vec[1]),
                    (combn(partition_vec[3],2)+sum(partition_vec[1:2])))
idx_null <- combn_null[1,]+num_partition*combn_null[2,]
combn_mat <- combn(num_partition,2)
idx_all <- combn_mat[1,]+num_partition*combn_mat[2,]

# # true postive rate (our method)
# idx <- which(idx_all %in% idx_null)
#
# hyp_tpr_null <- sapply(pval_vec_null, function(x){ #trial
#   pval <- p.adjust(x, method = "bonferroni")
#   alpha_vec <- seq(0, 1, length.out = 21)
#
#   sapply(alpha_vec, function(alpha){ #alpha
#     z <- which(pval >= alpha)
#     length(which(idx %in% z))/length(idx_null)
#   })
# })
# hyp_tpr_null <- apply(hyp_tpr_null, 1, mean)
#
# hyp_fpr_null <- sapply(pval_vec_null, function(x){ #trial
#   pval <- p.adjust(x, method = "bonferroni")
#   alpha_vec <- seq(0, 1, length.out = 21)
#
#   sapply(alpha_vec, function(alpha){ #alpha
#     z <- which(pval >= alpha)
#     length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
#   })
# })
# hyp_fpr_null <- apply(hyp_fpr_null, 1, mean)
#
# ##

hyp_tpr_alt <- sapply(pval_vec_alt, function(x){ #trial
  pval <- p.adjust(x, method = "bonferroni")
  alpha_vec <- seq(0, 1, length.out = 21)

  sapply(alpha_vec, function(alpha){ #alpha
    z <- which(pval >= alpha)
    length(which(idx %in% z))/length(idx_null)
  })
})
hyp_tpr_alt <- apply(hyp_tpr_alt, 1, mean)

hyp_fpr_alt <- sapply(pval_vec_alt, function(x){ #trial
  pval <- p.adjust(x, method = "bonferroni")
  alpha_vec <- seq(0, 1, length.out = 21)

  sapply(alpha_vec, function(alpha){ #alpha
    z <- which(pval >= alpha)
    length(which(!z %in% idx))/(length(idx_all) - length(idx_null))
  })
})
hyp_fpr_alt <- apply(hyp_fpr_alt, 1, mean)

par(mfrow = c(1,1))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
#lines(c(1,hyp_fpr_null,0), c(1,hyp_tpr_null,0), col = "red")
lines(c(1,hyp_fpr_alt,0), c(1,hyp_tpr_alt,0), col = "black")
lines(c(0,1), c(0,1), col = "red")



