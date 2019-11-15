color_palatte <- c(rgb(245, 234, 204, maxColorValue = 255), #yellow
                   rgb(189, 57, 60, maxColorValue = 255)) #red

# plot the graph
n <- length(dat_list)
g_selected <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g_selected <- igraph::add_edges(g_selected, edges = combn_mat[,stepdown_res$null_idx])

tmp <- rep(1, igraph::vcount(g_selected)); tmp[idx_our] <- 2
igraph::V(g_selected)$color <- color_palatte[tmp]
igraph::V(g_selected)$size <- c(5,10)[tmp]

# first construct 3 sets of nodes: first is our selected idx, the second is all the other nodes in the
## the giant component, and the third is all the remaining nodes
n <- igraph::vcount(g_selected)
idx1 <- intersect(idx_our, grep("PFC\\.[3-5]", names(dat_list)))
idx2 <- sort(idx_our)
idx2 <- setdiff(idx2, idx1)
tmp <-  igraph::components(g_selected)
idx3 <- sort(setdiff(which(tmp$membership == 1), c(idx1,idx2)))
adj_tmp <- as.matrix(igraph::as_adjacency_matrix(g_selected))
adj_tmp <- adj_tmp[c(idx1, idx2, idx3), c(idx1, idx2, idx3)]
diag(adj_tmp) <- 1

.rotate = function(a) { t(a[nrow(a):1,]) } 

png("../figures/figure_9.png", height = 1000, width = 2000, units = "px", res = 300)
par(mar = c(0,0,3,0), mfrow = c(1,2))
set.seed(10)
igraph::plot.igraph(g_selected, vertex.label = NA, main = "Full graph")

par(mar = c(3,3,3,0.5))
# next plot the adjacency matrix
image(.rotate(adj_tmp), asp = T, col = color_palatte, breaks = c(-.5,.5,1.5), xaxt = "n", yaxt = "n",
      xlab = "", ylab = "", main = "Adjacency matrix (subgraph)", axes = F)
title(ylab = "Index locations", mgp = c(1,1,0))
title(xlab = "Index locations", mgp = c(1,1,0))

# put in dashed lines
x_width <- length(idx_our)/nrow(adj_tmp)
y_height <- 1 - x_width
lines(rep(x_width, 2), c(1,1-x_width), lwd = 2, lty = 2)
lines(c(0,x_width), rep(y_height, 2), lwd = 2, lty = 2)
graphics.off()

############################################

png("../figures/figure_10b.png", height = 1300, width = 1150, units = "px", res = 300)
par(mar = c(4,4,4,1))
goodness_our <- goodness_list[[2]]
plot(sort(goodness_our), seq(0,1,length.out = length(goodness_our)), asp = T, xlab = "Theoretical quantiles",
     pch = 16,
     ylab = "Observed quantiles", main = "QQ-plot based on data\nfrom selected partitions")
lines(c(0,1),c(0,1), col = "red", lty = 2)
points(sort(goodness_our), seq(0,1,length.out = length(goodness_our)), pch = 16)
graphics.off()

##############################################

png("../figures/figure_2.png", height = 1300, width = 2300, units = "px", res = 300)
par(mar = c(4,4,4,1), mfrow = c(1,2))
plot(sort(goodness_pfc35), seq(0,1,length.out = length(goodness_pfc35)), asp = T, xlab = "Theoretical quantiles",
     pch = 16,
     ylab = "Observed quantiles", main = "QQ-plot based on partitions\nfrom only Window 1B")
lines(c(0,1),c(0,1), col = "red", lty = 2)
points(sort(goodness_pfc35), seq(0,1,length.out = length(goodness_pfc35)), pch = 16)

plot(sort(goodness_all), seq(0,1,length.out = length(goodness_all)), asp = T, xlab = "Theoretical quantiles",
     pch = 16,
     ylab = "Observed quantiles", main = "QQ-plot based on all partitions")
lines(c(0,1),c(0,1), col = "red", lty = 2)
points(sort(goodness_all), seq(0,1,length.out = length(goodness_all)), pch = 16)
graphics.off()

#########################################

color_palatte <- c(rgb(245, 234, 204, maxColorValue = 255), #yellow
                   rgb(189, 57, 60, maxColorValue = 255), #red
                   rgb(149,220,144, max = 255)) #green
library(igraph)

validated_genes <- covarianceSelection::validated_genes$Gene
g_our <- igraph::graph_from_adjacency_matrix(adj_our)

tmp <- rep(1, igraph::vcount(g_our))
tmp[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 3
tmp[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_our)$color <- color_palatte[tmp]
igraph::V(g_our)$size <- c(1,5,5)[tmp]
group_vec <- rep(NA, igraph::vcount(g_our))
group_vec[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 1
group_vec[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_our)$group <- group_vec
tmp <- igraph::components(g_our)
g_our <- igraph::induced_subgraph(g_our, v = which(tmp$membership == 1))
g_our <- igraph::as.undirected(g_our)
g_our <- igraph::minimum.spanning.tree(g_our)

png("../figures/tmp_1.png", height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(0,0,4,0))
set.seed(10)
igraph::plot.igraph(g_our, vertex.label = NA, main = "Our graph", edge.width = 0.5)
graphics.off()

idx1 <- which(igraph::V(g_our)$group == 1)
idx2 <- which(igraph::V(g_our)$group == 2)
dist_our <- igraph::distances(g_our, v = idx1, to = idx2)
set.seed(10)
tmp_graph <- igraph::erdos.renyi.game(igraph::vcount(g_our), igraph::ecount(g_our), type = "gnm")
dist_tmp <- igraph::distances(tmp_graph, v = 1:length(idx1), to = (length(idx1)+1):(length(idx1)+length(idx2)))
for(k in 1:10){
   tmp1 <- apply(dist_our, 1, function(x){sort(x, decreasing = F)[k]})
   #print(table(tmp1))
   print(mean(tmp1))
   #tmp2 <- mean(apply(dist_tmp, 1, function(x){sort(x, decreasing = F)[k]}))
   #print(tmp1/tmp2)
}

ratio_vec_our <- sapply(idx1, function(x){
   tmp <- as.numeric(igraph::neighborhood(g_our, order = 2, nodes = x)[[1]])
   length(intersect(tmp, idx2))/length(tmp)
})
mean(ratio_vec_our)

## 

g_pfc35 <- igraph::graph_from_adjacency_matrix(adj_pfc35)

tmp <- rep(1, igraph::vcount(g_pfc35))
tmp[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 3
tmp[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_pfc35)$color <- color_palatte[tmp]
igraph::V(g_pfc35)$size <- c(1,5,5)[tmp]
group_vec <- rep(NA, igraph::vcount(g_pfc35))
group_vec[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 1
group_vec[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_pfc35)$group <- group_vec
tmp <- igraph::components(g_pfc35)
g_pfc35 <- igraph::induced_subgraph(g_pfc35, v = which(tmp$membership == 1))
g_pfc35 <- igraph::as.undirected(g_pfc35)
g_pfc35 <- igraph::minimum.spanning.tree(g_pfc35)

png("../figures/tmp_2.png", height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(0,0,4,0))
set.seed(10)
igraph::plot.igraph(g_pfc35, vertex.label = NA, main = "PFC graph", edge.width = 0.5)
graphics.off()

idx1 <- which(igraph::V(g_pfc35)$group == 1)
idx2 <- which(igraph::V(g_pfc35)$group == 2)
dist_pfc35 <- igraph::distances(g_pfc35, v = idx1, to = idx2)
set.seed(10)
tmp_graph <- igraph::erdos.renyi.game(igraph::vcount(dist_pfc35), igraph::ecount(dist_pfc35), type = "gnm")
dist_tmp <- igraph::distances(tmp_graph, v = 1:length(idx1), to = (length(idx1)+1):(length(idx1)+length(idx2)))
for(k in 1:10){
   tmp1 <- apply(dist_pfc35, 1, function(x){sort(x, decreasing = F)[k]})
   print(mean(tmp1))
   #tmp2 <- mean(apply(dist_tmp, 1, function(x){sort(x, decreasing = F)[k]}))
   #print(tmp1/tmp2)
}

ratio_vec_pfc35 <- sapply(idx1, function(x){
   tmp <- as.numeric(igraph::neighborhood(g_pfc35, order = 2, nodes = x)[[1]])
   length(intersect(tmp, idx2))/length(tmp)
})
mean(ratio_vec_pfc35)

######################################

library(igraph)

k_vec <- 1:10
res_mat <- matrix(NA, ncol = 2, nrow = length(k_vec))
validated_genes <- covarianceSelection::validated_genes$Gene

g_our <- igraph::graph_from_adjacency_matrix(adj_our)
tmp <- rep(1, igraph::vcount(g_our))
group_vec <- rep(NA, igraph::vcount(g_our))
group_vec[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 1
group_vec[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_our)$group <- group_vec
tmp <- igraph::components(g_our)
g_our <- igraph::induced_subgraph(g_our, v = which(tmp$membership == 1))
g_our <- igraph::as.undirected(g_our)
g_our <- igraph::minimum.spanning.tree(g_our)

idx1 <- which(igraph::V(g_our)$group == 1)
idx2 <- which(igraph::V(g_our)$group == 2)
dist_tmp <- igraph::distances(tmp_graph, v = 1:length(idx1), to = (length(idx1)+1):(length(idx1)+length(idx2)))
for(k in k_vec){
   tmp1 <- apply(dist_our, 1, function(x){sort(x, decreasing = F)[k]})
   res_mat[k,1] <- mean(tmp1)
}

g_pfc35 <- igraph::graph_from_adjacency_matrix(adj_pfc35)
group_vec <- rep(NA, igraph::vcount(g_pfc35))
group_vec[which(colnames(dat_list[[1]]) %in% validated_genes)] <- 1
group_vec[which(colnames(dat_list[[1]]) %in% genes_nodawn)] <- 2
igraph::V(g_pfc35)$group <- group_vec
tmp <- igraph::components(g_pfc35)
g_pfc35 <- igraph::induced_subgraph(g_pfc35, v = which(tmp$membership == 1))
g_pfc35 <- igraph::as.undirected(g_pfc35)
g_pfc35 <- igraph::minimum.spanning.tree(g_pfc35)

idx1 <- which(igraph::V(g_pfc35)$group == 1)
idx2 <- which(igraph::V(g_pfc35)$group == 2)
dist_tmp <- igraph::distances(tmp_graph, v = 1:length(idx1), to = (length(idx1)+1):(length(idx1)+length(idx2)))
for(k in k_vec){
   tmp1 <- apply(dist_pfc35, 1, function(x){sort(x, decreasing = F)[k]})
   res_mat[k,2] <- mean(tmp1)
}


plot(NA, xlim = range(k_vec), ylim = range(res_mat))
points(res_mat[,1])
points(res_mat[,2], col = "red")

#########################
res_mat <- matrix(NA, ncol = 2, nrow = 100)

for(k in 1:100){
   our_mat <- eigen_our$vectors[,c(1:k, (ncol(eigen_our$vectors)-k+1):ncol(eigen_our$vectors))]
   idx1 <- which(colnames(dat_list[[1]]) %in% validated_genes)
   idx2 <- which(colnames(dat_list[[1]]) %in% genes_nodawn)
   idx1 <- setdiff(idx1, idx2)
   
   dist_mat <- sapply(idx1, function(x){
      sapply(idx2, function(y){
         .l2norm(our_mat[x,] - our_mat[y,])
      })
   })
   val1 <- mean(apply(dist_mat, 1, min))
   
   ########
   
   pfc_mat <- eigen_pfc35$vectors[,c(1:k, (ncol(eigen_pfc35$vectors)-k+1):ncol(eigen_pfc35$vectors))]
   idx1 <- which(colnames(dat_list[[1]]) %in% validated_genes)
   idx2 <- which(colnames(dat_list[[1]]) %in% genes_nodawn)
   idx1 <- setdiff(idx1, idx2)
   
   dist_mat <- sapply(idx1, function(x){
      sapply(idx2, function(y){
         .l2norm(pfc_mat[x,] - our_mat[y,])
      })
   })
   val2 <- mean(apply(dist_mat, 1, min))
   
   res_mat[k,] <- c(val1, val2)
}

plot(NA, xlim = c(1,100), ylim = range(res_mat))
points(res_mat[,1])
points(res_mat[,2], col = "red")

