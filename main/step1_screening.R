p_thres_screening <- 0.1 
cor_thres_screening <- 0.75

###

if(verbose) print(paste0(Sys.time(), "Start of step 1: Screening"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 13964

screening_res <- covarianceSelection::screen(dat_pfc35, pv = tada$pval.TADA, p_thres = p_thres_screening, 
                                             cor_thres = cor_thres_screening)

# 1605 primary, 1460 secondary, total of 3065

# apply the new gene list
for(i in 1:length(dat_list)){
  dat_list[[i]] <- dat_list[[i]][,c(screening_res$primary, screening_res$secondary)]
}
tada <- tada[c(screening_res$primary, screening_res$secondary),]
                   
if(verbose) print(paste0("Dimension of dat_list is: ", unique(sapply(dat_list, ncol)), collapse = ", "))

rm(list = c("selected_idx", "dat_pfc35", "i"))

save.image(file = paste0(save_filepath, "/step1_screening.RData"))