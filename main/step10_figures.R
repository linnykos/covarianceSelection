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
# set.seed(10)
# idx3 <- sample(idx3)
adj_tmp <- as.matrix(igraph::as_adjacency_matrix(g_selected))
# idx3 <- idx3[order(rowSums(adj_tmp[idx3,]), decreasing = T)]
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
