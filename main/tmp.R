load("../results/step7_results.RData")

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::clique_selection(g_sub, threshold = 0.90)[[1]]]
zz <- clique_selection(g_selected, threshold = 0.9, target_idx = core_set, verbose = T, max_length = 10000)

covarianceSelection::binning(names(dat_list)[zz[[1]]])

idx_our <- zz[[1]]
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)

res <- covarianceSelection::graphicalModel(dat_our, primary_idx = 1:length(screening_res$primary), lambda = 0.075)
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

validated_genes <- covarianceSelection::validated_genes$Gene
num_our <- length(intersect(genes_our, validated_genes))
length(validated_genes)
num_our
