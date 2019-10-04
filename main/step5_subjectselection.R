trials <- 100

if(verbose) print(paste0(Sys.time(), "Start of step 5: Subject selection"))

save(trials, file = paste0(save_filepath, "/test.RData"))
stepdown_obj <- covarianceSelection::stepdown_path(dat_list, trials = trials, cores = ncores, verbose = verbose,
                                          iterations = 7, prob = 1-1e-5,
                                          file = paste0(save_filepath, "/step5_subjectselection_tmp2.RData"))
save.image(file = paste0(save_filepath, "/step5_subjectselection_tmp.RData"))
alternative_test_vec <- covarianceSelection::stepdown(dat_list, prob = 0.99, only_test_stat = T)
save.image(file = paste0(save_filepath, "/step5_subjectselection_tmp.RData"))
stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  covarianceSelection::stepdown_choose(stepdown_obj, alpha = alpha, return_pvalue = T)
})

save.image(file = paste0(save_filepath, "/step5_subjectselection.RData"))

######
# sapply(stepdown_res, function(x){
#   length(x$null_idx)
# })