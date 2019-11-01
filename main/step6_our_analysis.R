gamma_threshold <- 0.95

if(verbose) print(paste0(Sys.time(), "Start of step 6: Our data analysis"))

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::clique_selection(g_sub, threshold = gamma_threshold)[[1]]]
idx_our <- covarianceSelection::clique_selection(g_selected, threshold = gamma_threshold, target_idx = core_set)
idx_our <- idx_our[[1]]

## covarianceSelection::binning(names(dat_list)[idx_our]); covarianceSelection::binning(names(dat_list))
## igraph::ecount(igraph::induced_subgraph(g_selected, idx_our))/(length(idx_our)*(length(idx_our)-1)/2)
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

# res <- covarianceSelection::graphicalModel_range(dat_our, 1:length(screening_res$primary), lambda_min = 0.01, lambda_max = 0.35,  
#                                                  lambda_length = 30, verbose = T) 
# save.image(file = paste0(save_filepath, "/step6_ourdata_analysis", filepath_suffix, ".RData"))
# 
# scale_vec_our <- sapply(res, function(x){covarianceSelection::compute_scale_free(as.matrix(x$adj_mat))})
# edges_vec_our <- sapply(res, function(x){sum(as.matrix(x$adj_mat))/2})
# # idx <- which.max(scale_vec_our)
# idx <- 26
# adj_our <- as.matrix(res[[idx]]$adj_mat)
# stopifnot(all(dim(adj_our) == nrow(tada)))

res <- covarianceSelection::graphicalModel(dat_our, primary_idx = 1:length(screening_res$primary), lambda =  seq(0.05, 0.1, length.out = 15)[5])
adj_our <- as.matrix(res$adj_mat)

# run the HMRF
set.seed(10)
seedindex <- rep( 0, ncol(adj_our))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
set.seed(10)
hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres) 
report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= fdr_cutoff)]))

adj_our <- Matrix::Matrix(adj_our, sparse = T)

rm(list = c("dat_our", "seedindex", "res", "combn_mat", "n", "g_selected"))

save.image(file = paste0(save_filepath, "/step6_ourdata_analysis", filepath_suffix, ".RData"))