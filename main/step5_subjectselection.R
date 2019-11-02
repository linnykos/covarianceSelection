trials <- 100
prob_val <- 1-1e-5
stepdown_alpha <- 0.1

if(verbose) print(paste0(Sys.time(), "Start of step 5: Subject selection"))

# save(trials, file = paste0(save_filepath, "/test.RData"))
# stepdown_obj <- covarianceSelection::stepdown_path(dat_list, trials = trials, cores = ncores, verbose = verbose,
#                                           iterations = 7, prob = prob_val,
#                                           file = paste0(save_filepath, "/step5_subjectselection_tmp2", filepath_suffix, ".RData"))
# save.image(file = paste0(save_filepath, "/step5_subjectselection_tmp", filepath_suffix, ".RData"))
# alternative_test_vec <- covarianceSelection::stepdown(dat_list, prob = 0.99, only_test_stat = T)
# save.image(file = paste0(save_filepath, "/step5_subjectselection_tmp", filepath_suffix, ".RData"))
# stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
#   covarianceSelection::stepdown_choose(stepdown_obj, alpha = alpha, return_pvalue = T)
# })

set.seed(10)
stepdown_res <- covarianceSelection::stepdown(dat_list, trials = trials, alpah = stepdown_alpha, return_pvalue = T,
                                              prob = prob_val, verbose = T)

save.image(file = paste0(save_filepath, "/step5_subjectselection", filepath_suffix, ".RData"))

######
# sapply(stepdown_res, function(x){
#   length(x$null_idx)
# })