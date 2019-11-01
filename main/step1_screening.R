p_thres_screening <- 0.01 
p_primary <- 0.1
num_genes <- 3500

###

if(verbose) print(paste0(Sys.time(), "Start of step 1: Screening"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 13964

screening_res <- covarianceSelection::screen(dat_pfc35, pv = tada$pval.TADA, p_thres = p_thres_screening, 
                                             num_genes = num_genes)

# 265 primary, 3235 secondary, total of 3500

# reorder which genes are primary and which are secondary
all_idx <- sort(unique(c(screening_res$primary, screening_res$secondary)))
screening_res$primary <- all_idx[which(tada$pval.TADA[all_idx] < p_primary)]
screening_res$secondary <- setdiff(all_idx, screening_res$primary)

# apply the new gene list
for(i in 1:length(dat_list)){
  dat_list[[i]] <- dat_list[[i]][,c(screening_res$primary, screening_res$secondary)]
}
tada <- tada[c(screening_res$primary, screening_res$secondary),]
                
if(verbose) print(paste0("Dimension of dat_list is: ", unique(sapply(dat_list, ncol)), collapse = ", "))

rm(list = c("selected_idx", "dat_pfc35", "i"))

save.image(file = paste0(save_filepath, "/step1_screening", filepath_suffix, ".RData"))