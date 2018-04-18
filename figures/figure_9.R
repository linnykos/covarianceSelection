rm(list=ls())
library(igraph)

load("../results/step2_res.RData")
load("../results/step5_res.RData")
.l2norm <- function(x){sqrt(sum(x^2))}

extra_nodes <- 10

set.seed(10)
adj <- as.matrix(igraph::as_adj(g))
diag(adj) <- 1
idx <- sample(selected_idx)

deg_vec <- igraph::degree(g)
idx2 <- c(1:len)[-idx]
deg_vec2 <- deg_vec[idx2]
idx2 <- idx2[order(deg_vec2, decreasing = F)[floor(length(deg_vec)/2.25):(floor(length(deg_vec)/2.25)+extra_nodes)]]

adj_small <- adj[c(idx, idx2), c(idx, idx2)]

.clockwise90 = function(a) { t(a[nrow(a):1,]) }

#colors
pale <- rgb(247, 234, 200, max = 255)
red <- rgb(205,40,54,maxColorValue=255)

pdf("../figures/figure_9.pdf", height = 3, width = 6.25)
par(mar = c(3,3,2,1), mfrow = c(1,2))
image(.clockwise90(adj_small), asp = T, xlab = "Index locations", ylab = "Index locations",
      main = "Adjacency matrix (Subgraph)",
      xaxt = "n", bty = "n", yaxt = "n", col = c(pale, red), mgp = c(1,0,0))

spacing <- 2*1/(2*(length(idx)+extra_nodes)-2)
len_base <- length(idx)
lines(x = c(-.5*spacing, (len_base-1.5)*spacing), y = rep(1+(1.5-len_base)*spacing, 2), lwd = 2, lty = 2)
lines(x = rep((len_base-1.5)*spacing, 2), y = c(1+spacing*.5, 1+(1.5-len_base)*spacing), lwd = 2, lty = 2)

#################
#par(mar = c(0, 0, 4, 0))
node_col <- rep(pale, igraph::vcount(g))
node_size <- rep(4, igraph::vcount(g))
node_col[selected_idx] <- red
node_size[selected_idx] <- 6

set.seed(10)
l <- igraph::layout_nicely(g)

#manipulate graph size
center <- apply(l[selected_idx,], 2, mean)
comp <- igraph::components(g)
idx <- which(comp$membership == which.max(comp$csize))
dist_mat1 <- apply(l[idx,], 1, function(x){
  .l2norm(x - center)
})
dist_mat2 <- apply(l[-idx,], 1, function(x){
  .l2norm(x - center)
})

ratio <- max(dist_mat1)/min(dist_mat2)*1.1

idx2 <- c(1:nrow(l))[-idx]
for(i in idx2){
  d <- l[i,] - center
  d <- d*ratio
  l[i,] <- center + d
}

igraph::plot.igraph(g, layout = l, vertex.size = node_size, vertex.label = NA,
                    vertex.color = node_col, main = "Full graph")

graphics.off()
