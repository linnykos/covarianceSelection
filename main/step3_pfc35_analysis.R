fdr_cutoff <- 0.01
pthres <- 0.05

#####
set.seed(10)
if(verbose) print(paste0(Sys.time(), "Start of step 3: Naive analysis"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
# res <- covarianceSelection::graphicalModel_range(dat_pfc35, 1:length(screening_res$primary), 
#                                                  lambda_min = 0.01, lambda_max = 0.35, 
#                                                  lambda_length = 30, verbose = T) 
# save.image(file = paste0(save_filepath, "/step3_pfc35_analysis", filepath_suffix, ".RData"))
# 
# scale_vec_pfc35 <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
# edges_vec_pfc35 <- sapply(res, function(x){sum(as.matrix(x$adj_mat))/2})
# # idx <- which.max(scale_vec_pfc35)
# idx <- 27
# adj_pfc35 <- as.matrix(res[[idx]]$adj_mat)

res <- covarianceSelection::graphicalModel(dat_pfc35, primary_idx = 1:length(screening_res$primary), lambda = seq(0.05, 0.1, length.out = 15)[1])
adj_pfc35 <- as.matrix(res$adj_mat)
stopifnot(all(dim(adj_pfc35) == nrow(tada)))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
set.seed(10)
hmrf_pfc35 <- covarianceSelection::hmrf(tada$pval.TADA, adj_pfc35, seedindex, pthres = pthres)
report_pfc35 <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35$post, tada$pval.TADA, hmrf_pfc35$Iupdate)
genes_pfc35 <- sort(as.character(report_pfc35$Gene[which(report_pfc35$FDR <= fdr_cutoff)]))

adj_pfc35 <- Matrix::Matrix(adj_pfc35, sparse = T)

rm(list = c("dat_pfc35", "seedindex", "idx", "res", "selected_idx"))

save.image(file = paste0(save_filepath, "/step3_pfc35_analysis", filepath_suffix, ".RData"))

########################

# plot(exp(seq(log(0.01), log(0.35), length.out = length(scale_vec_pfc35))), scale_vec_pfc35)
