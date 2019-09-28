rm(list=ls())
set.seed(10)
n <- 30
combn_mat <- utils::combn(n, 2)
edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
g <- igraph::graph.empty(n = n, directed = F)
g <- igraph::add_edges(g, edges = edge_mat)
# res <- chen_2010(g)

#########
threshold = 0.95

g <- igraph::as.undirected(g)
g <- igraph::simplify(g)

n <- igraph::vcount(g)
if(igraph::ecount(g)/choose(n,2) >= threshold) return(as.numeric(igraph::V(g)$name))
igraph::V(g)$name <- 1:n
c_matrix <- .form_c_matrix(g)

q <- dequer::deque()
dequer::push(q, .chen_object(g, c_matrix))
max_size <- 0
max_node_set <- numeric(0)

while(length(q) > 0){
  obj <- dequer::pop(q)
  n_internal <- igraph::vcount(obj$g)
  if(n_internal <= max_size) next()
  
  node_set_internal <- as.character(igraph::V(obj$g)$name)
  den <- .chen_check_density(g, node_set_internal)
  if(den >= threshold){
    max_size <- n_internal
    max_node_set <- node_set_internal
  }
  
  if(igraph::ecount(obj$g) == 0 | nrow(obj$c_matrix) == 0) next()
  
  res <- .chen_separate(obj$g, obj$c_matrix, threshold = threshold)
  dequer::push(q, .chen_object(res$g1, res$c_matrix1))
  dequer::push(q, .chen_object(res$g2, res$c_matrix2))
}

#####################
g <- obj$g
c_matrix <- obj$c_matrix
check = T

n <- igraph::vcount(g)
node_set <- as.numeric(igraph::V(g)$name)
stopifnot(length(node_set) == n)
if(check)  stopifnot(all(c_matrix[,1] %in% node_set), all(c_matrix[,2] %in% node_set))

idx <- 1
comp_res <- igraph::components(g)

while(comp_res$no == 1){
  edge_id <- igraph::get.edge.ids(g, which(node_set %in% c_matrix[idx,1:2]))
  if(edge_id != 0){
    g <- igraph::delete.edges(g, edge_id)
    comp_res <- igraph::components(g)
  } 
  
  idx <- idx + 1
}

idx1 <- which(comp_res$membership == 1)
idx2 <- c(1:n)[-idx1]
stopifnot(length(idx1) > 0 & length(idx2) > 0)

g1 <- igraph::induced_subgraph(g, idx1)
g2 <- igraph::induced_subgraph(g, idx2)

c_matrix1 <- c_matrix[intersect(which(c_matrix[,1] %in% node_set[idx1]), 
                                which(c_matrix[,2] %in% node_set[idx1])),,drop = F]
c_matrix2 <- c_matrix[intersect(which(c_matrix[,1] %in% node_set[idx2]), 
                                which(c_matrix[,2] %in% node_set[idx2])),,drop = F]
