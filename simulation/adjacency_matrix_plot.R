#produces figures of the adjacency matrix

rm(list=ls())
load("../results/nonparanormal.RData")

idx <- 3
n <- sum(paramMat[1,1:3])

# grab the desired edges
index <- res[[idx]][[1]]$indices_list[[20]] 

# form the graph
g <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(n, 2)
edges <- combn_mat[,index]
g <- igraph::add_edges(g, edges = edges)

# form partition
tmp <- covarianceSelection::clique_selection(g, threshold = 0.95)
partition <- tmp[[1]]

#compute the adjacency matrix
adj <- matrix(0, n, n)
diag(adj) <- 1
for(i in 1:ncol(edges)){
  adj[edges[1,i], edges[2,i]] <- 1
  adj[edges[2,i], edges[1,i]] <- 1
}

#rearrange so the first 5 rows are the ones we start with, shuffle the rest
set.seed(10)
base_idx <- c(1:3, 16, 21)
idx <- sample(c(1:n)[-base_idx])
adj2 <- adj[c(base_idx, idx), c(base_idx, idx)]
adj3 <- adj[c(partition, c(1:25)[-partition]), c(partition, c(1:25)[-partition])]

.clockwise90 = function(a) { t(a[nrow(a):1,]) }

#colors
pale <- rgb(247, 234, 200, max = 255)
red <- rgb(205,40,54,maxColorValue=255)

png("../figures/figure_3.png", height = 1200, width = 2400, res = 300, units = "px")
#plot the first image
par(mfrow = c(1,2), mar = c(3,3,2,1))
image(.clockwise90(adj2), asp = T, xlab = "Index locations", ylab = "Index locations",
      main = "Adjacency matrix",
      xaxt = "n", bty = "n", yaxt = "n", col = c(pale, red), mgp = c(1,0,0))

spacing <- 2*1/(2*25-2)
len_base <- length(base_idx)
# lines(x = c(-.5*spacing, (len_base-.5)*spacing), y = rep(1+(.5-len_base)*spacing, 2), lwd = 2, lty = 2)
# lines(x = rep((len_base-.5)*spacing, 2), y = c(1+spacing*.5, 1+(.5-len_base)*spacing), lwd = 2, lty = 2)

#plot the second image

image(.clockwise90(adj3), asp = T, xlab = "Index locations", ylab = "Index locations",
      main = "Adjacency matrix (reordered)",col = c(pale, red), mgp = c(1,0,0),
      xaxt = "n", bty = "n", yaxt = "n")

len <- length(partition)
lines(x = c(-.5*spacing, (len-.5)*spacing), y = rep(1+(.5-len)*spacing, 2), lwd = 2, lty = 2)
lines(x = rep((len-.5)*spacing, 2), y = c(1+spacing*.5, 1+(.5-len)*spacing), lwd = 2, lty = 2)
graphics.off()
