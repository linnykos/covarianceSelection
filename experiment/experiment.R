rm(list=ls())
ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection.RData")

combn_mat <- utils::combn(length(dat_list), 2)
combn_mat <- combn_mat[,sample(1:ncol(combn_mat), 1000)]
pval_vec <- sapply(1:ncol(combn_mat), function(i){
  x <- combn_mat[,i]
  res <- covarianceSelection::cai_test(dat_list[[x[1]]], dat_list[[x[2]]], cores = ncores, prob = 0.99)
  print(paste0("Indicies ", i, " : value ", round(res, 3)))
  res
})

save.image("/raid6/Kevin/covarianceSelection/results/tmp.RData")