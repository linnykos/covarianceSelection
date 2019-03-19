rm(list=ls())
load("../results/clique_simulation_small.RData")

index <- res[[2]][[1]]$indices_list[[13]]
n <- 14
edges <- utils::combn(n, 2)

g <- igraph::graph.empty(n = n, directed = F)
g <- igraph::add_edges(g, edges = edges[, index])
