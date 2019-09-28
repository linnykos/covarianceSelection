rm(list=ls())
load("../results/step4_ouranalysis_tmp.RData")

res <- tsourakakis_2013(g)
###########
threshold = 0.95
iter_max = round(igraph::vcount(g)/2)

g <- igraph::as.undirected(g)
g <- igraph::simplify(g)
igraph::V(g)$name <- 1:n

initial_idx <- as.character(.tsourakakis_initialize(g))
node_set <- sort(c(as.character(igraph::neighbors(g, v = initial_idx)), initial_idx))
iter <- 1
# print(node_set)

while(iter <= iter_max){
  while(TRUE){
    # print(node_set)
    den_org <- .tsourakakis_obj(g, threshold, node_set)
    node_candidate <- setdiff(as.character(igraph::V(g)$name), node_set)
    # print(node_candidate)
    # print(class(node_candidate))
    next_set <- .find_candidate(g, threshold, node_set, node_candidate, den_org)
    if(any(is.na(next_set))) break()
    node_set <- next_set
  }
  
  den_org <- .tsourakakis_obj(g, threshold, node_set)
  next_set <- .find_candidate(g, threshold, node_set, NA, den_org)
  if(any(is.na(next_set))) break()
  node_set <- next_set
  iter <- iter+1
  # print(node_set)
}