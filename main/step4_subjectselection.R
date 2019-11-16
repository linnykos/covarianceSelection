trials <- 100
prob_val <- 1-1e-5
stepdown_alpha <- 0.1

if(verbose) print(paste0(Sys.time(), "Start of step 5: Subject selection"))

set.seed(10)
stepdown_res <- covarianceSelection::stepdown(dat_list, trials = trials, alpha = stepdown_alpha, return_pvalue = F,
                                              prob = prob_val, verbose = T, cores = ncores)

save.image(file = paste0(save_filepath, "/step4_subjectselection", filepath_suffix, ".RData"))
