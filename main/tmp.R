rm(list=ls())
load("/raid6/Kevin/covarianceSelection/results/step7_results.RData")

ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)

res_pfc35 <- covarianceSelection::graphicalModel_range(dat_pfc35, 1:length(screening_res$primary), 
                                                 lambda_min = 0.05, lambda_max = 0.1, lambda_length = 15, verbose = F)

validated_genes <- covarianceSelection::validated_genes$Gene

for(i in 1:length(res_pfc35)){
  adj_pfc35 <- as.matrix(res_pfc35[[i]]$adj_mat)
  
  # run the HMRF
  set.seed(10)
  seedindex <- rep(0, ncol(adj_pfc35))
  seedindex[which(tada$dn.LoF >= 3)] <- 1
  
  set.seed(10)
  hmrf_pfc35 <- covarianceSelection::hmrf(tada$pval.TADA, adj_pfc35, seedindex, pthres = pthres)
  report_pfc35 <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35$post, tada$pval.TADA, hmrf_pfc35$Iupdate)
  # genes_pfc35 <- sort(as.character(report_pfc35$Gene[which(report_pfc35$FDR <= fdr_cutoff)]))
  genes_pfc35 <- sort(as.character(report_pfc35$Gene[order(report_pfc35$FDR)[1:200]]))
  
  num_pfc35 <- length(intersect(genes_pfc35, validated_genes)) 
  num_edges <- sum(adj_pfc35)/2
  scale_free <- covarianceSelection::compute_scale_free(adj_pfc35)
  
  print(paste0("Level: ", i, " // Edges: ", num_edges, " // Scale-free: ", round(scale_free,2), 
               " // Total: ", length(genes_pfc35), " // Num: ", num_pfc35))
}

#####################################

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::clique_selection(g_sub, threshold = 0.95)[[1]]]
idx_our <- covarianceSelection::clique_selection(g_selected, threshold = 0.95, target_idx = core_set, verbose = T, max_length = 5000)
idx_our <- idx_our[[1]]

covarianceSelection::binning(names(dat_list)[idx_our]); covarianceSelection::binning(names(dat_list))
igraph::ecount(igraph::induced_subgraph(g_selected, idx_our))/(length(idx_our)*(length(idx_our)-1)/2)

dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

res_our <- covarianceSelection::graphicalModel_range(dat_our, 1:length(screening_res$primary), 
                                                       lambda_min = 0.05, lambda_max = 0.1, lambda_length = 15, verbose = F)


for(i in 1:length(res_our)){
  adj_our <- as.matrix(res_our[[i]]$adj_mat)
  
  # run the HMRF
  set.seed(10)
  seedindex <- rep(0, ncol(adj_our))
  seedindex[which(tada$dn.LoF >= 3)] <- 1
  
  set.seed(10)
  hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres)
  report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
  #genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= fdr_cutoff)]))
  genes_our <- sort(as.character(report_our$Gene[order(report_our$FDR)[1:200]]))
  
  num_our <- length(intersect(genes_our, validated_genes)) 
  num_edges <- sum(adj_our)/2
  scale_free <- covarianceSelection::compute_scale_free(adj_our)
  
  print(paste0("Level: ", i, " // Edges: ", num_edges, " // Scale-free: ", round(scale_free,2), 
               " // Total: ", length(genes_our), " // Num: ", num_our))
}

#########################

# robustness

gamma_vec <- seq(0.85, 0.97, length.out = 13)
gene_list <- vector("list", length = length(gamma_vec))

for(i in 1:length(gamma_vec)){
  n <- length(dat_list)
  g_selected <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(length(dat_list), 2)
  g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res$null_idx])
  
  # construct the core set
  selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
  g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
  core_set <- selected_idx[covarianceSelection::clique_selection(g_sub, threshold = gamma_vec[i])[[1]]]
  idx_our <- covarianceSelection::clique_selection(g_selected, threshold = gamma_vec[i], 
                                                   target_idx = core_set, verbose = T, max_length = 5000)
  idx_our <- idx_our[[1]]
  
  print(covarianceSelection::binning(names(dat_list)[idx_our]))
  # igraph::ecount(igraph::induced_subgraph(g_selected, idx_our))/(length(idx_our)*(length(idx_our)-1)/2)
  
  dat_our <- do.call(rbind, dat_list[idx_our])
  dat_our <- scale(dat_our, scale = F)
  
  res <- covarianceSelection::graphicalModel(dat_our, primary_idx = 1:length(screening_res$primary), lambda = seq(0.05, 0.1, length.out = 15)[5])
  adj_our <- as.matrix(res$adj_mat)
  
  # run the HMRF
  set.seed(10)
  seedindex <- rep(0, ncol(adj_our))
  seedindex[which(tada$dn.LoF >= 3)] <- 1
  
  set.seed(10)
  hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, adj_our, seedindex, pthres = pthres)
  report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
  genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= fdr_cutoff)]))
  # genes_our <- sort(as.character(report_our$Gene[order(report_our$FDR)[1:200]]))
  gene_list[[i]] <- genes_our
  
  num_our <- length(intersect(genes_our, validated_genes)) 
  num_edges <- sum(adj_our)/2
  scale_free <- covarianceSelection::compute_scale_free(adj_our)
  
  print(paste0("Gamma: ", gamma_vec[i], " // Edges: ", num_edges, " // Scale-free: ", round(scale_free,2), 
               " // Total: ", length(genes_our), " // Num: ", num_our))
}

gene_intersect <- intersect(gene_list[[1]], gene_list[[2]])
for(i in 3:length(gene_list)){
  gene_intersect <- intersect(gene_intersect, gene_list[[3]])
}
length(intersect(gene_intersect, validated_genes)) 
  
  