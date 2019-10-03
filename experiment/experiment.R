rm(list=ls())
load("../results/step3_pfc35_analysis.RData")
sapply(res, function(x){compute_scale_free(as.matrix(x$adj_mat))})

adj_pfc35 <- as.matrix(res[[15]]$adj_mat)

# run the HMRF
set.seed(10)
seedindex <- rep(0, nrow(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_pfc35 <- covarianceSelection::hmrf(tada$pval.TADA, adj_pfc35, seedindex, pthres = pthres)
report_pfc35 <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35$post, tada$pval.TADA, hmrf_pfc35$Iupdate)
cutoff <- sort(report_pfc35$FDR, decreasing = FALSE)[num_target]
genes_pfc35 <- sort(as.character(report_pfc35$Gene[which(report_pfc35$FDR <= cutoff)]))

hmrf_pfc35$c # seems to be around -2??
table(report_pfc35$indicator)
quantile(report_pfc35$FDR)

validated_genes <- read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))
validated_genes <- covarianceSelection::symbol_synonyms(validated_genes, verbose = T)

length(intersect(genes_pfc35, validated_genes))

###################
load("../results/step4_alldata_analysis.RData")
sapply(res, function(x){compute_scale_free(as.matrix(x$adj_mat))})

adj_all <- as.matrix(res[[14]]$adj_mat)

set.seed(10)
seedindex <- rep(0, ncol(adj_all))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_all <- covarianceSelection::hmrf(tada$pval.TADA, adj_all, seedindex, pthres = pthres) 
report_all <- covarianceSelection::report_results(tada$Gene, 1-hmrf_all$post, tada$pval.TADA, hmrf_all$Iupdate)
cutoff <- sort(report_all$FDR, decreasing = FALSE)[num_target]
genes_all <- sort(as.character(report_all$Gene[which(report_all$FDR <= cutoff)]))

hmrf_all$c # seems to be around 2??
table(report_all$indicator)
quantile(report_all$FDR)

length(intersect(genes_all, validated_genes)) 
