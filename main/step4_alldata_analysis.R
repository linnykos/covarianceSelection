if(verbose) print(paste0(Sys.time(), "Start of step 3: All data analysis"))

dat_all <- do.call(rbind, dat_list) # 107 x 3438
dat_all <- scale(dat_all, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection:::graphicalModel_range(dat_all, 1:length(screening_res$primary), lambda_min = 0.01, lambda_max = 0.35, verbose = T) 
save.image(file = paste0(save_filepath, "/step4_alldata_analysis.RData"))

scale_vec_all <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
idx <- which.max(scale_vec_all)
adj_all <- as.matrix(res[[idx]]$adj_mat)
stopifnot(all(dim(adj_all) == nrow(tada)))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_all))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_all <- covarianceSelection::hmrf(tada$pval.TADA, adj_all, seedindex, pthres = pthres) 
report_all <- covarianceSelection::report_results(tada$Gene, 1-hmrf_all$post, tada$pval.TADA, hmrf_all$Iupdate)
cutoff <- sort(report_all$FDR, decreasing = FALSE)[num_target]
genes_all <- sort(as.character(report_all$Gene[which(report_all$FDR <= cutoff)]))

rm(list = c("dat_all", "seedindex", "cutoff", "res", "scale_idx", "idx"))

save.image(file = paste0(save_filepath, "/step4_alldata_analysis.RData"))