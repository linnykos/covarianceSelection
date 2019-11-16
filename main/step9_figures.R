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

validated_idx <- which(colnames(dat_list[[1]]) %in% validated_genes)
nodawn_idx <- which(colnames(dat_list[[1]]) %in% genes_nodawn)
validated_idx <- setdiff(validated_idx, nodawn_idx)
k_vec <- 1:10

res_mat1 <- sapply(k_vec, function(k){
   val1 <- covarianceSelection::compute_mst_distance(as.matrix(adj_our), validated_idx, nodawn_idx, k)
   val2 <- covarianceSelection::compute_mst_distance(as.matrix(adj_pfc35), validated_idx, nodawn_idx, k)
   
   c(val1, val2)
})

k_vec <- 1:10
res_mat2 <- sapply(k_vec, function(k){
   val1 <- covarianceSelection::compute_graph_root_distance(eigen_our$vectors, validated_idx, nodawn_idx, k)
   val2 <- covarianceSelection::compute_graph_root_distance(eigen_pfc35$vectors, validated_idx, nodawn_idx, k)
   
   c(val1, val2)
})

(res_mat2[2,]-res_mat2[1,])/res_mat2[1,]

png("../figures/appendix_9.png", height = 1300, width = 2300, units = "px", res = 300)
par(mar = c(4,4,4,1), mfrow = c(1,2))
plot(NA, xlim = range(k_vec), ylim = range(res_mat1), main = "Distance from Satterstrom genes to\nDe Rubeis genes: MST",
     xlab = "Number of closeby De Rubeis genes", ylab = "Distance in MST", cex.main = 1)
lines(k_vec, res_mat1[2,], lwd = 2, col = color_palatte[1])
points(k_vec, res_mat1[2,], pch = 21, bg = color_palatte[1])
lines(k_vec, res_mat1[1,], lwd = 2, col = color_palatte[2])
points(k_vec, res_mat1[1,], pch = 16, col = color_palatte[2])
legend("topleft", c("Using Window 1B graph", "Using COBS graph"),
       bty="n", fill=color_palatte, cex = 0.75)

plot(NA, xlim = range(k_vec), ylim = range(res_mat2), main = "Distance from Satterstrom genes to\nDe Rubeis genes: graph root embedding",
     xlab = "Number of closeby De Rubeis genes", ylab = "Distance in graph root embedding", cex.main = 1)
lines(k_vec, res_mat2[2,], lwd = 2, col = color_palatte[1])
points(k_vec, res_mat2[2,], pch = 21, bg = color_palatte[1])
lines(k_vec, res_mat2[1,], lwd = 2, col = color_palatte[2])
points(k_vec, res_mat2[1,], pch = 16, col = color_palatte[2])

legend("topleft", c("Using Window 1B graph", "Using COBS graph"),
       bty="n", fill=color_palatte, cex = 0.75)
graphics.off()

median((res_mat1[2,] - res_mat1[1,])/res_mat1[2,])
median((res_mat2[2,] - res_mat2[1,])/res_mat2[2,])
