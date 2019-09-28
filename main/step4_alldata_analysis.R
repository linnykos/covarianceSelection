if(verbose) print(paste0(Sys.time(), "Start of step 3: All data analysis"))

dat_all <- do.call(rbind, dat_list) # 107 x 3438
dat_all <- scale(dat_all)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
adj_all <- covarianceSelection::graphicalModel(dat_all, lambda = "lambda.1se", verbose = T) 
stopifnot(all(dim(adj_all) == nrow(tada)))

save.image(file = paste0(save_filepath, "/step4_alldata_analysis.RData"))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_all))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_all <- covarianceSelection::hmrf(tada$pval.TADA, adj_all, seedindex, pthres = pthres) #??? WARNING???
report_all <- covarianceSelection::report_results(tada$Gene, 1-hmrf_all$post, tada$pval.TADA, hmrf_all$Iupdate)
cutoff <- sort(report_all$FDR, decreasing = FALSE)[num_target]
genes_all <- sort(as.character(report_all$Gene[which(report_all$FDR <= cutoff)]))

rm(list = c("dat_all", "seedindex", "cutoff"))

save.image(file = paste0(save_filepath, "/step4_alldata_analysis.RData"))