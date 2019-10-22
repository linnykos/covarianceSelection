selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection::graphicalModel_range(dat_pfc35, 1:length(screening_res$primary), 
                                                 lambda_min = 0.01, lambda_max = 0.35, 
                                                 lambda_length = 30, verbose = T) 

validated_genes <- covarianceSelection::validated_genes$Gene

vec <- sapply(20:length(res), function(x){
  print(x)
  
  # run the HMRF
  set.seed(10)
  seedindex <- rep(0, ncol(adj_pfc35))
  seedindex[which(tada$dn.LoF >= 3)] <- 1
  
  set.seed(10)
  hmrf_pfc35 <- covarianceSelection::hmrf(tada$pval.TADA, as.matrix(res[[x]]$adj_mat), seedindex, pthres = pthres)
  report_pfc35 <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35$post, tada$pval.TADA, hmrf_pfc35$Iupdate)
  genes_pfc35 <- sort(as.character(report_pfc35$Gene[which(report_pfc35$FDR <= fdr_cutoff)]))
  
  length(intersect(genes_pfc35, validated_genes))
})

################################

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::tsourakakis_2013(g_sub)]

idx_our <- covarianceSelection::chen_2010(g_selected, core_set = core_set)
## covarianceSelection::binning(names(dat_list)[idx_our]); covarianceSelection::binning(names(dat_list))
## igraph::ecount(igraph::induced_subgraph(g_selected, idx_our))/(length(idx_our)*(length(idx_our)-1)/2)
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

res <- covarianceSelection::graphicalModel_range(dat_our, 1:length(screening_res$primary), lambda_min = 0.01, lambda_max = 0.35,  
                                                 lambda_length = 30, verbose = T) 

vec3 <- sapply(20:length(res), function(x){
  print(x)
  
  # run the HMRF
  set.seed(10)
  seedindex <- rep(0, ncol(adj_pfc35))
  seedindex[which(tada$dn.LoF >= 3)] <- 1
  
  set.seed(10)
  hmrf_our <- covarianceSelection::hmrf(tada$pval.TADA, as.matrix(res[[x]]$adj_mat), seedindex, pthres = pthres)
  report_our <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our$post, tada$pval.TADA, hmrf_our$Iupdate)
  genes_our <- sort(as.character(report_our$Gene[which(report_our$FDR <= fdr_cutoff)]))
  
  length(intersect(genes_our, validated_genes))
})
