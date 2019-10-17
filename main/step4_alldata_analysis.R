if(verbose) print(paste0(Sys.time(), "Start of step 3: All data analysis"))

dat_all <- do.call(rbind, dat_list) # 107 x 3438
dat_all <- scale(dat_all, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection::graphicalModel_range(dat_all, 1:length(screening_res$primary), lambda_min = 0.01, lambda_max = 0.35,  
                                                 lambda_length = 30, verbose = T) 
save.image(file = paste0(save_filepath, "/step4_alldata_analysis", filepath_suffix, ".RData"))

scale_vec_all <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
edges_vec_all <- sapply(res, function(x){sum(as.matrix(x$adj_mat))/2})
# idx <- which.max(scale_vec_all)
idx <- 21
adj_all <- as.matrix(res[[idx]]$adj_mat)
stopifnot(all(dim(adj_all) == nrow(tada)))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_all))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
set.seed(10)
hmrf_all <- covarianceSelection::hmrf(tada$pval.TADA, adj_all, seedindex, pthres = pthres) 
report_all <- covarianceSelection::report_results(tada$Gene, 1-hmrf_all$post, tada$pval.TADA, hmrf_all$Iupdate)
genes_all <- sort(as.character(report_all$Gene[which(report_all$FDR <= fdr_cutoff)]))

adj_all <- Matrix::Matrix(adj_all, sparse = T)

rm(list = c("dat_all", "seedindex", "res", "idx"))

save.image(file = paste0(save_filepath, "/step4_alldata_analysis", filepath_suffix, ".RData"))

########################

# plot(exp(seq(log(0.01), log(0.35), length.out = length(scale_vec_all))), scale_vec_all)

