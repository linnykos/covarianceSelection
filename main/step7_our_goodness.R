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
goodness_our1 <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 250, prob = 1-1e-5/2)

save.image(file = paste0(save_filepath, "/step7_our_goodness", filepath_suffix, ".RData"))

goodness_our2 <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 250, prob = 1-1e-5)

save.image(file = paste0(save_filepath, "/step7_our_goodness", filepath_suffix, ".RData"))

goodness_our3 <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 250, prob = 1-1e-4)

save.image(file = paste0(save_filepath, "/step7_our_goodness", filepath_suffix, ".RData"))

# hist(goodness_our, col = "gray", breaks = 20)
# plot(sort(goodness_our), seq(0,1,length.out = length(goodness_our)), asp = T); lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = "red")
