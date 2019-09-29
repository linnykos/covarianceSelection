if(verbose) print(paste0(Sys.time(), "Start of step 6: Our data analysis"))

idx_our <- covarianceSelection:::tsourakakis_2013(g)
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our)
adj_our <- covarianceSelection::graphicalModel(dat_our, lambda = "lambda.1se", verbose = T) 

save.image(file = paste0(save_filepath, "/step6_ourdata_analysis.RData"))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_our))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres) #??? WARNING???
report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
cutoff <- sort(report_our$FDR, decreasing = FALSE)[num_target]
genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= cutoff)]))

save.image(file = paste0(save_filepath, "/step6_ourdata_analysis.RData"))