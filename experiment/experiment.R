rm(list=ls())
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection.RData")
prob = 0.9999

dat_list <- lapply(dat_list, scale, center = T, scale = F)
diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
len <- length(dat_list)
combn_mat <- utils::combn(len, 2)

num_list <- lapply(dat_list, function(x){covarianceSelection:::.compute_sigma(x, diag_idx)})
denom_list <- covarianceSelection:::.compute_all_denom(dat_list, num_list, diag_idx)

ncores <- 20
doMC::registerDoMC(cores = ncores)

func <- function(x){
  print(x)
  covarianceSelection:::.compute_covStat(num_list[[combn_mat[1,x]]], num_list[[combn_mat[2,x]]],
                   denom_list[[combn_mat[1,x]]], denom_list[[combn_mat[2,x]]],
                   squared = T, prob = prob)
}

t_vec <- foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), func(i))

save(t_vec, "/raid6/Kevin/covarianceSelection/results/step3_subjectselection_updated.RData")

################################