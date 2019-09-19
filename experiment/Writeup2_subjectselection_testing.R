# in this test, we see if reducing the rank of the covariance helps at all
rm(list=ls())

ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
#install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)

verbose <- T
save_filepath <- "/raid6/Kevin/covarianceSelection/results"

load("../results/step1_screening.RData")

##########################

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- dat_list[selected_idx]

# rescale them, and keep only the first 3 principal components
k <- 3
dat_pfc35 <- lapply(dat_pfc35, function(x){
  scale(x, center = T, scale = T)
})

# dat_pfc35_2 <- lapply(dat_pfc35, function(x){
#   svd_res <- svd(x)
#   res <- svd_res$u[,1:k] %*% diag(svd_res$d[1:k]) %*% t(svd_res$v[,1:k])
#   scale(res, center = T, scale = T)
# })

dat_list <- dat_pfc35
dat_list <- lapply(dat_list, scale, center = T, scale = F)
diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
len <- length(dat_list)
combn_mat <- utils::combn(len, 2)

num_list <- lapply(dat_list, function(x){.compute_sigma(x, diag_idx)})
denom_list <- .compute_all_denom(dat_list, num_list, diag_idx)
num_x <- num_list[[1]]
num_y <- num_list[[3]]
denom_x <- denom_list[[1]]
denom_y <- denom_list[[3]]
zz <- (num_x - num_y)^2/(denom_x + denom_y)

############################

ncores <- 20
trials <- 200

save(trials, file = paste0(save_filepath, "/test.RData"))
stepdown_obj <- covarianceSelection::stepdown_path(dat_pfc35, trials = trials, cores = ncores, verbose = verbose,
                                                   iterations = 7, file = paste0(save_filepath, "/step3_subjectselection_tmp2.RData"),
                                                   prob = 0.9999)
save.image(file = paste0(save_filepath, "/step3_subjectselection_tmp.RData"))
stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  covarianceSelection::stepdown_choose(stepdown_obj, alpha = alpha, return_pvalue = T)
})

save.image(file = paste0(save_filepath, "/step3_subjectselection_experiment_3.RData"))
