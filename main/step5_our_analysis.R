gamma_threshold <- 0.95

if(verbose) print(paste0(Sys.time(), "Start of step 6: Our data analysis"))

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::clique_selection(g_sub, threshold = gamma_threshold)[[1]]]
idx_our <- covarianceSelection::clique_selection(g_selected, threshold = gamma_threshold, target_idx = core_set)
idx_our <- idx_our[[1]]

dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

res <- covarianceSelection::graphicalModel(dat_our, primary_idx = 1:length(screening_res$primary), 
                                           lambda = seq(0.05, 0.1, length.out = 15)[5])
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

rm(list = c("dat_our", "seedindex", "res", "combn_mat", "n", "g_selected", "g_sub", "selected_idx", "core_set"))

save.image(file = paste0(save_filepath, "/step5_ourdata_analysis", filepath_suffix, ".RData"))



