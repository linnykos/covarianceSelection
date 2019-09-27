rm(list=ls())
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection_updated2.RData")

doMC::registerDoMC(cores = 20)

idx <- covarianceSelection:::tsourakakis_2014_approximate(g)
dat_list <- lapply(dat_list[idx], scale)
dat_our <- do.call(rbind, dat_list)
adj_our <- covarianceSelection::graphicalModel(dat_our, lambda = "lambda.1se", verbose = T) 

save.image(file = paste0(save_filepath, "/step4_ouranalysis_tmp.RData"))

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(adj_our))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres) #??? WARNING???
report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
cutoff <- sort(report_our$FDR, decreasing = FALSE)[num_target]
genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= cutoff)]))

save.image(file = paste0(save_filepath, "/step4_ouranalysis_tmp.RData"))