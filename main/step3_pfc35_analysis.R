pthres <- 0.05
num_target <- 200

#####
set.seed(10)
if(verbose) print(paste0(Sys.time(), "Start of step 2: Naive analysis"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3438
dat_pfc35 <- scale(dat_pfc35, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection:::graphicalModel_range(dat_pfc35, screening_res$primary, lambda_min = 0.01, lambda_max = 0.35, verbose = T) 

save.image(file = paste0(save_filepath, "/step3_pfc35_analysis.RData"))
scale_idx <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
idx <- which.max(scale_idx)
adj_pfc35 <- as.matrix(res[[idx]]$adj_mat)
stopifnot(all(dim(adj_pfc35) == nrow(tada)))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_pfc35 <- covarianceSelection::hmrf(tada$pval.TADA, adj_pfc35, seedindex, pthres = pthres)
report_pfc35 <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35$post, tada$pval.TADA, hmrf_pfc35$Iupdate)
cutoff <- sort(report_pfc35$FDR, decreasing = FALSE)[num_target]
genes_pfc35 <- sort(as.character(report_pfc35$Gene[which(report_pfc35$FDR <= cutoff)]))

rm(list = c("dat_pfc35", "seedindex", "cutoff", "scale_idx", "idx", "res", "selected_idx"))

save.image(file = paste0(save_filepath, "/step3_pfc35_analysis.RData"))