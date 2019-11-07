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

# goodness_our <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 100, prob = 1,
#                                                            verbose = T)
# print(paste0("Finished our"))

# dat_our <- do.call(rbind, dat_list[idx_our])
# dat_our <- scale(dat_our, scale = F)
# 
prob_vec <- c(1, 1-1e-5, 1-1e-4)
goodness_list <- vector("list", length = length(prob_vec))

for(i in 1:length(prob_vec)){
  set.seed(10)
  goodness_list[[i]] <- covarianceSelection::goodness_of_fit(dat_list[idx_our], permutations = 250, trials = 100, prob = prob_vec[i],
                                                             verbose = T)
  print(paste0("Goodness ", i, " done"))
  save.image(file = paste0(save_filepath, "/step7_our_goodness", filepath_suffix, ".RData"))
}

save.image(file = paste0(save_filepath, "/step7_our_goodness", filepath_suffix, ".RData"))

# hist(goodness_our, col = "gray", breaks = 20)
# plot(sort(goodness_list[[1]]), seq(0,1,length.out = length(goodness_list[[1]])), asp = T); lines(c(0,1),c(0,1), col = "red", lty = 2)
