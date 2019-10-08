if(verbose) print(paste0(Sys.time(), "Start of step 6: Our data analysis"))

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

idx_our <- covarianceSelection:::tsourakakis_2013(g_selected)
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

res <- covarianceSelection:::graphicalModel_range(dat_our, 1:length(screening_res$primary), lambda_min = 0.01, lambda_max = 0.35, verbose = T) 
save.image(file = paste0(save_filepath, "/step6_ourdata_analysis.RData"))

scale_vec_our <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
idx <- which.max(scale_vec_our)
adj_our <- as.matrix(res[[idx]]$adj_mat)
stopifnot(all(dim(adj_our) == nrow(tada)))

# run the HMRF
set.seed(10)
seedindex <- rep( 0, ncol(adj_our))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres) #??? WARNING???
report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
cutoff <- sort(report_our$FDR, decreasing = FALSE)[num_target]
genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= cutoff)]))

rm(list = c("dat_our", "seedindex", "cutoff", "res", "combn_mat", "n", "g_selected"))

save.image(file = paste0(save_filepath, "/step6_ourdata_analysis.RData"))