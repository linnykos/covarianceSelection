trials <- 100

len <- length(dat_list)

stepdown_obj <- covarianceSelection::stepdown_path(dat_list, trials = trials, cores = ncores, verbose = verbose,
                                          iterations = 7, file = paste0(save_filepath, "/step3_subjectselection_tmp2.RData"))
save.image(file = paste0(save_filepath, "/step3_subjectselection_tmp.RData"))
stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  covarianceSelection::stepdown_choose(stepdown_obj, alpha = alpha, return_pvalue = T)
})

save.image(file = paste0(save_filepath, "/step3_subjectselection.RData"))