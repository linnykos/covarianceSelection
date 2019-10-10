adj_pfc35 = as.matrix(adj_pfc35)
idx = which(colSums(adj_pfc35) == 0)
adj = adj_pfc35[-idx, -idx]
g = igraph::graph_from_adjacency_matrix(adj)

png("../figures/test.png", height = 1500, width = 1500, res = 300, units = "px")
plot(g)
graphics.off()
