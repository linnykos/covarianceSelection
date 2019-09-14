trials <- 200

len <- length(dat_list)

stepdown_obj <- covarianceSelection::stepdown_path(dat_list, trials = trials, cores = ncores, verbose = verbose,
                                          iterations = 10)
stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  covarianceSelection::stepdown_choose(obj, alpha = alpha, return_pvalue = T)
})

save.image(file = paste0(save_filepath, "/step3_subjectselection.RData"))