pthres_screening <- 0.1 
num_genes <- 3500

###

if(verbose) print(paste0(Sys.time(), "Start of step 1: Screening"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 13964

screening_res <- covarianceSelection::screen(dat_pfc35, pv = tada$pval.TADA, pthres = pthres_screening, num_genes = num_genes)
# 1605 primary, 1895 secondary

all_genes <- sort(unique(c(screening_res$primary, screening_res$secondary))) # 3500 genes

# apply the new gene list
for(i in 1:length(dat_list)){
  dat_list[[i]] <- dat_list[[i]][,all_genes]
}
tada <- tada[all_genes,]
                   
if(verbose) print(paste0("Dimension of dat_list is: ", unique(sapply(dat_list, ncol)), collapse = ", "))

rm(list = c("selected_idx", "dat_pfc35", "v_seq", "res_pca", "screening_res", "i"))

save.image(file = paste0(save_filepath, "/step1_screening.RData"))