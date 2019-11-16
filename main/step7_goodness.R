if(verbose) print(paste0(Sys.time(), "Start of step 7: Computing goodness of fit"))
prob_goodness <- seq(1-1e-5, 1-1e-4, length.out = 5)[2]

goodness_our <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 100, verbose = T, prob = prob_goodness)

if(verbose) print(paste0(Sys.time(), "Finished ours"))
save.image(file = paste0(save_filepath, "/step7_goodness", filepath_suffix, ".RData"))

##

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
set.seed(10)
goodness_pfc35 <- covarianceSelection::goodness_of_fit(dat_list[selected_idx], permutations = 250, trials = 100, 
                                                       verbose = T, prob = prob_goodness)

if(verbose) print(paste0(Sys.time(), "Finished PFC35"))
save.image(file = paste0(save_filepath, "/step7_goodness", filepath_suffix, ".RData"))

##

set.seed(10)
prob_goodness <- 1-1e-6
goodness_all <- covarianceSelection::goodness_of_fit(dat_list, permutations = 250, trials = 100, prob = prob_goodness, 
                                                     verbose = T)

if(verbose) print(paste0(Sys.time(), "Finished all"))
rm(list = c("prob_goodness"))

save.image(file = paste0(save_filepath, "/step7_goodness", filepath_suffix, ".RData"))
