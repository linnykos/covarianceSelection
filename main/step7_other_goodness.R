selected_idx <- grep("PFC\\.[3-5]", names(dat_list))

set.seed(10)
goodness_pfc35 <- covarianceSelection::goodness_of_fit(dat_list[selected_idx], permutations = 250, trials = 100, 
                                                       verbose = T)

save.image(file = paste0(save_filepath, "/step8_other_goodness", filepath_suffix, ".RData"))
print(paste0("Finished pfc35"))

# dat_all <- do.call(rbind, dat_list) 
# dat_all <- scale(dat_all, scale = F)

set.seed(10)
goodness_val <- 1-1e-6
goodness_all <- covarianceSelection::goodness_of_fit(dat_list, permutations = 250, trials = 100, prob = goodness_val, 
                                                     verbose = T)
print(paste0("Finished all"))

save.image(file = paste0(save_filepath, "/step8_other_goodness", filepath_suffix, ".RData"))

# hist(goodness_pfc35, col = "gray", breaks = 20); hist(goodness_all, col = "gray", breaks = 20)
# plot(sort(goodness_pfc35), seq(0,1,length.out = length(goodness_pfc35)), asp = T); lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = "red")
# plot(sort(goodness_all), seq(0,1,length.out = length(goodness_all)), asp = T); lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = "red")


