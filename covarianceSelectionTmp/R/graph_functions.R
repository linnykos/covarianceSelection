#' Compute the scale-free statistic
#'
#' @param adj adjacency matrix
#'
#' @return numeric
#' @export
compute_scale_free <- function(adj){
  stopifnot(is.matrix(adj), is.numeric(adj), nrow(adj) == ncol(adj))

  adj[which(adj != 0, arr.ind = TRUE)] <- 1

  tmp <- apply(adj,1,sum)
  x <- table(tmp)+1
  y <- as.numeric(names(x))+1

  stats::cor(log(x),log(y))^2
}

#' List the index-pairs of non-zeros in an adjacency matrix
#'
#' Only lists elements in the upper-triange, excluding diagonal
#'
#' @param adj matrix
#' @param tol small number
#'
#' @return matrix with 2 columns
#' @export
enumerate_edges <- function(adj, tol = 1e-6){
  stopifnot(sum(is.na(adj))==0, ncol(adj) == nrow(adj))
  stopifnot(is.matrix(adj) | class(adj) == "dgCMatrix")

  diag(adj) <- 0
  adj[lower.tri(adj)] <- 0
  which(abs(adj) > tol, arr.ind = TRUE)
}

#' Convert adjancency matrix to graph
#'
#' @param adj adjacency matrix
#' @param tol tolerance
#'
#' @return \code{igraph} object
#' @export
adj_to_graph <- function(adj, tol = 1e-6){
  n <- ncol(adj)
  idx <- which(abs(adj) >= tol, arr.ind = T)
  g <- igraph::graph.empty(n = n, directed = F)
  igraph::add_edges(g, edges = t(idx))
}

#' Count neighbors of a graph
#' 
#' Count how many neighbors each vertex in \code{idx} has.
#'
#' @param g igraph object
#' @param idx vector of indicies
#' @param order what degree neighbor
#'
#' @return a vector of numerics
#' @export
count_neighbors <- function(g, idx, order = 1) {
  stopifnot(class(g) == "igraph")
  stopifnot(all(idx %% 1 == 0) & max(idx) <= igraph::vcount(g))

  neigh <- igraph::ego(g, order = order)
  neigh_vec <- numeric(igraph::vcount(g))

  for(i in 1:length(neigh)) {
    #make sure it doesn't count itself
    neigh_vec[i] = length(which(idx %in% setdiff(neigh[[i]],i)))
  }

  neigh_vec
}

