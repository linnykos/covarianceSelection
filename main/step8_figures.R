color_palatte <- c(rgb(245, 234, 204, maxColorValue = 255), #yellow
                   rgb(189, 57, 60, maxColorValue = 255)) #red

# plot the graph
n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

tmp <- rep(1, igraph::vcount(g_selected)); tmp[idx_our] <- 2
igraph::V(g_selected)$color <- color_palatte[tmp]

set.seed(10)
igraph::plot.igraph(g_selected, vertex.label = NA, vertex.size = 10)

# zoom in on the big component
tmp <-  igraph::components(g_selected)
tmp_g <- igraph::induced_subgraph(g_selected, which(tmp$membership == 1))
set.seed(10)
igraph::plot.igraph(tmp_g, vertex.label = NA, vertex.size = 10)

# try another strategy, plot the subgraph containing only neighbors of pfc35 nodes
idx_pfc35 <- grep("PFC\\.[3-5]", names(dat_list))
idx_neigh <- unique(c(unlist(sapply(idx_pfc35, function(x){
  igraph::neighbors(g_selected, x)
})), idx_pfc35))
tmp_g <- igraph::induced_subgraph(g_selected, idx_neigh)
set.seed(10)
igraph::plot.igraph(tmp_g, vertex.label = NA, vertex.size = 10)

#################

# plot the adjacency matrix

# first construct 3 sets of nodes: first is our selected idx, the second is all the other nodes in the
## the giant component, and the third is all the remaining nodes
n <- igraph::vcount(g_selected)
idx1 <- intersect(idx_our, grep("PFC\\.[3-5]", names(dat_list)))
idx2 <- sort(idx_our)
tmp <-  igraph::components(g_selected)
idx3 <- sort(setdiff(which(tmp$membership == 1), c(idx1,idx2)))
idx4 <- sort(setdiff(1:n, c(idx1,idx2,idx3)))
adj_tmp <- as.matrix(igraph::as_adjacency_matrix(g_selected))
adj_tmp <- adj_tmp[c(idx1, idx2, idx3, idx4), c(idx1, idx2, idx3, idx4)]

.rotate = function(a) { t(a[nrow(a):1,]) } 
image(.rotate(adj_tmp), asp = T, col = color_palatte, breaks = c(-.5,.5,1.5))

# put in dashed lines
x_width <- length(idx_our)/n
y_height <- 1 - x_width
lines(rep(x_width, 2), c(0,1), lwd = 2, lty = 2)
lines(c(0,1), rep(y_height, 2), lwd = 2, lty = 2)

##########################

plot(report_pfc35$FDR, report_our$FDR, asp = T)

###########################

# plot p-value vs the number of risk neighbors
adj_pfc35 <- as.matrix(adj_pfc35)
d <- ncol(adj_pfc35)
diag(adj_pfc35) <- 0
idx_gene <- which(colnames(dat_list[[1]]) %in% genes_pfc35)
vec_pfc35 <- sapply(1:d, function(i){
  if(i %% floor(ncol(adj_pfc35)/10) == 0) cat('*')
  neigh_vec <- which(adj_pfc35[,i] == 1)
  length(intersect(neigh_vec, idx_gene))
})
col_vec <- rep(1, d); col_vec[idx_gene] <- 2
plot(tada$pval.TADA, jitter(vec_pfc35), col = color_palatte[col_vec], pch = 16)

###################
# 
# # trying something janky
# idx <- which(colSums(adj_pfc35) != 0)
# tada_zz <- tada[idx,]
# adj_zz <- adj_pfc35[idx, idx]
# 
# set.seed(10)
# seedindex <- rep(0, ncol(adj_pfc35))
# seedindex[which(tada$dn.LoF >= 3)] <- 1
# seedindex <- seedindex[idx]
# 
# set.seed(10)
# hmrf_pfc35_zz <- covarianceSelection::hmrf(tada_zz$pval.TADA, adj_zz, seedindex, pthres = pthres)
# plot(tada_zz$pval.TADA, 1-hmrf_pfc35_zz$post, asp = T)

idx <- which(colnames(dat_list[[1]]) %in% validated_genes)
colSums(adj_pfc35)[idx]
colSums(adj_our)[idx]

###############################
load("../results/step6_ourdata_analysis_3531.RData")

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3065
dat_pfc35 <- scale(dat_pfc35, scale = F)
dim(dat_pfc35)
res <- covarianceSelection::graphicalModel(dat_pfc35, primary_idx = 1:length(screening_res$primary), lambda = 0.075)
adj_pfc35_zz <- as.matrix(res$adj_mat)
covarianceSelection::compute_scale_free(adj_pfc35_zz) #0.84
sum(adj_pfc35_zz/2) #6884
table(colSums(adj_pfc35_zz))
set.seed(10)
seedindex <- rep(0, ncol(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1
hmrf_pfc35_zz <- hmrf(tada$pval.TADA, adj_pfc35_zz, seedindex, pthres = pthres)
hmrf_pfc35_zz$c
report_pfc35_zz <- covarianceSelection::report_results(tada$Gene, 1-hmrf_pfc35_zz$post, tada$pval.TADA, hmrf_pfc35_zz$Iupdate)
genes_pfc35_zz <- sort(as.character(report_pfc35_zz$Gene[which(report_pfc35_zz$FDR <= fdr_cutoff)]))
plot(report_pfc35_zz$FDR, tada$pval.TADA, asp = T)
length(intersect(genes_pfc35_zz, validated_genes)) 
length(genes_pfc35_zz) #34/185 
which(report_pfc35_zz$FDR <= fdr_cutoff)

idx <- which(colnames(dat_list[[1]]) %in% validated_genes)
colSums(adj_pfc35_zz)[idx]

d <- ncol(adj_pfc35_zz)
diag(adj_pfc35_zz) <- 0
idx_gene <- which(colnames(dat_list[[1]]) %in% genes_pfc35_zz)
vec_pfc35 <- sapply(1:d, function(i){
  if(i %% floor(ncol(adj_pfc35_zz)/10) == 0) cat('*')
  neigh_vec <- which(adj_pfc35_zz[,i] == 1)
  length(intersect(neigh_vec, idx_gene))
})
col_vec <- rep(1, d); col_vec[idx_gene] <- 2
plot(tada$pval.TADA, jitter(vec_pfc35), col = color_palatte[col_vec], pch = 16, xlim = c(0,0.4))
plot(tada$qvalue, jitter(vec_pfc35), col = color_palatte[col_vec], pch = 16)

### 

n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res[[3]]$null_idx])

# construct the core set
selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
g_sub <- igraph::induced_subgraph(g_selected, selected_idx)
core_set <- selected_idx[covarianceSelection::chen_2010(g_sub)]

idx_our <- covarianceSelection::chen_2010(g_selected, core_set = core_set)
covarianceSelection::binning(names(dat_list)[idx_our]); covarianceSelection::binning(names(dat_list))
igraph::ecount(igraph::induced_subgraph(g_selected, idx_our))/(length(idx_our)*(length(idx_our)-1)/2)
dat_our <- do.call(rbind, dat_list[idx_our])
dat_our <- scale(dat_our, scale = F)
dim(dat_our)
res <- covarianceSelection::graphicalModel(dat_our, primary_idx = 1:length(screening_res$primary), lambda = 0.075)
adj_our_zz <- as.matrix(res$adj_mat)
covarianceSelection::compute_scale_free(adj_our_zz) #0.88 
sum(adj_our_zz/2) #7612
set.seed(10)
seedindex <- rep(0, ncol(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1
hmrf_our_zz <- hmrf(tada$pval.TADA, adj_our_zz, seedindex, pthres = pthres)
hmrf_our_zz$c
report_our_zz <- covarianceSelection::report_results(tada$Gene, 1-hmrf_our_zz$post, tada$pval.TADA, hmrf_our_zz$Iupdate)
genes_our_zz <- sort(as.character(report_our_zz$Gene[which(report_our_zz$FDR <= fdr_cutoff)]))
plot(report_our_zz$FDR, tada$pval.TADA, asp = T)
length(intersect(genes_our_zz, validated_genes)) 
length(genes_our_zz) #36/157
which(report_our_zz$FDR <= fdr_cutoff)

idx <- which(colnames(dat_list[[1]]) %in% validated_genes)
colSums(adj_our_zz)[idx]

d <- ncol(adj_our_zz)
diag(adj_our_zz) <- 0
idx_gene <- which(colnames(dat_list[[1]]) %in% genes_our_zz)
vec_our <- sapply(1:d, function(i){
  if(i %% floor(ncol(adj_our_zz)/10) == 0) cat('*')
  neigh_vec <- which(adj_our_zz[,i] == 1)
  length(intersect(neigh_vec, idx_gene))
})
col_vec <- rep(1, d); col_vec[idx_gene] <- 2
plot(tada$pval.TADA, jitter(vec_our), col = color_palatte[col_vec], pch = 16, xlim = c(0,0.4))
plot(tada$qvalue, jitter(vec_our), col = color_palatte[col_vec], pch = 16)

###################


dat_all <- do.call(rbind, dat_list) # 107 x 3438
dat_all <- scale(dat_all, scale = F)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection::graphicalModel(dat_all, primary_idx = 1:length(screening_res$primary), lambda = 0.075)
adj_all_zz <- as.matrix(res$adj_mat)
covarianceSelection::compute_scale_free(adj_all_zz) #0.89
sum(adj_all_zz/2) #7424
set.seed(10)
seedindex <- rep(0, ncol(adj_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1
hmrf_all_zz <- hmrf(tada$pval.TADA, adj_all_zz, seedindex, pthres = pthres)
hmrf_all_zz$c
report_all_zz <- covarianceSelection::report_results(tada$Gene, 1-hmrf_all_zz$post, tada$pval.TADA, hmrf_all_zz$Iupdate)
genes_all_zz <- sort(as.character(report_all_zz$Gene[which(report_all_zz$FDR <= fdr_cutoff)]))
plot(report_all_zz$FDR, tada$pval.TADA, asp = T)
length(intersect(genes_all_zz, validated_genes)) 
length(genes_all_zz) #17/173
which(report_all_zz$FDR <= fdr_cutoff)

idx <- which(colnames(dat_list[[1]]) %in% validated_genes)
colSums(adj_all_zz)[idx]

d <- ncol(adj_all_zz)
diag(adj_all_zz) <- 0
idx_gene <- which(colnames(dat_list[[1]]) %in% genes_our_zz)
vec_our <- sapply(1:d, function(i){
  if(i %% floor(ncol(adj_all_zz)/10) == 0) cat('*')
  neigh_vec <- which(adj_all_zz[,i] == 1)
  length(intersect(neigh_vec, idx_gene))
})
col_vec <- rep(1, d); col_vec[idx_gene] <- 2
plot(tada$pval.TADA, jitter(vec_our), col = color_palatte[col_vec], pch = 16, xlim = c(0,0.4))
plot(tada$qvalue, jitter(vec_our), col = color_palatte[col_vec], pch = 16)
