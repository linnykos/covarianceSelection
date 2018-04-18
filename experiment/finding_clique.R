rm(list=ls())
library(longitudinalGM)
load("../longitudinalGM/tests/assets/clique_selection1.RData")

n <- 25
max_val <- 2^25

edges <- combn(25,2)[,lis[[1]]]
g <- igraph::graph.empty(n = n, directed = F)
g <- igraph::add_edges(g, edges = edges)
adj <- as.matrix(igraph::as_adjacency_matrix(g))

cores <- 14
breakpoints <- ceiling(seq(1, max_val, length.out = cores + 1))
dif <- max(diff(breakpoints))

func <- function(i){
  counter <- i
  max_clique <- NA

  while(counter <= min(i + dif, max_val)){
    if(i == 1 && counter %% floor(dif/10) == 0) print('*')

    vec <- binaryLogic::as.binary(counter)
    if(length(vec) <= n-1){
      vec <- c(rep(0, n-length(vec)), vec)
    }

  idx <- which(vec == 1)
    if(longitudinalGM:::.pass_threshold(adj[idx, idx, drop = F], 0.95)){
      if(length(idx) > length(max_clique)) max_clique <- idx
    }

    counter <- counter+1
  }

  max_clique
}

doMC::registerDoMC(cores = cores)
res <- foreach::"%dopar%"(foreach::foreach(i = 1:cores), func(breakpoints[i]))

save.image("finding_clique.RData")
