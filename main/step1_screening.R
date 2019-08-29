pthres_screening <- 0.1 

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 13962

primary_genes <- which(tada$pval.TADA <= pthres_screening) # 1605 genes

v_seq <- exp(seq(log(1), log(sqrt(ncol(dat_pfc35))), length.out = 10))
res_pca <- PMA::SPC(dat_pfc35, sumabsv = v_seq[5], K = 6, trace = F)
secondary_genes <- which(apply(res_pca$v, 1, function(x){any(abs(x) > 1e-3)})) # 2121 genes

all_genes <- sort(unique(c(primary_genes, secondary_genes))) # 3437 genes

# apply the new gene list
for(i in 1:length(dat_list)){
  dat_list[[i]] <- dat_list[[i]][,all_genes]
}
                   
rm(list = c("selected_idx", "dat_pfc35", "v_seq", "res_pca", 
            "primary_genes", "secondary_genes", "i"))

save.image(file = paste0(save_filepath, "/step1_res.RData"))