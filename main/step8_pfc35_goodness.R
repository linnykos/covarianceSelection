selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)

goodness_pfc35 <- covarianceSelection::goodness_of_fit(dat_pfc35, trials = 250, prob = 1-1e-5)

save.image(file = paste0(save_filepath, "/step8_pfc35_goodness", filepath_suffix, ".RData"))
