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


set.seed(10)
goodness_pfc35 <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 250, prob = prob_val)

save.image(file = paste0(save_filepath, "/step8_our_goodness", filepath_suffix, ".RData"))