rm(list=ls())
ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection.RData")

combn_mat <- utils::combn(length(dat_list), 2)
pval_vec <- apply(combn_mat, 2, function(x){
  print(x)
  covarianceSelection::cai_test(dat_list[[x[1]]], dat_list[[x[2]]], cores = ncores, prob = 0.99)
})

save.image("/raid6/Kevin/covarianceSelection/results/tmp.RData")