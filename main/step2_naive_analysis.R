if(verbose) print("Start of step 2: Naive analysis")

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3437

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
prec_mat <- covarianceSelection::graphicalModel(dat_pfc35, verbose = T) 

save.image(file = paste0(save_filepath, "/step2_naive_analysis.RData"))