selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)

goodness_pfc35 <- covarianceSelection::goodness_of_fit(dat_pfc35, permutations = 250, trials = 250, prob = 1)

save.image(file = paste0(save_filepath, "/step8_pfc35_goodness", filepath_suffix, ".RData"))

dat_all <- do.call(rbind, dat_list) # 107 x 3438
dat_all <- scale(dat_all, scale = F)

goodness_all <- covarianceSelection::goodness_of_fit(dat_all, permutations = 250, trials = 250, prob = 1)

save.image(file = paste0(save_filepath, "/step8_pfc35_goodness", filepath_suffix, ".RData"))

# hist(goodness_pfc35, col = "gray", breaks = 20)
# plot(sort(goodness_pfc35), seq(0,1,length.out = length(goodness_pfc35)), asp = T); lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = "red")