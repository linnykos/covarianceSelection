selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3437

prec_mat <- covarianceSelection::graphicalModel(dat_pfc35, verbose = T)
