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
  x <- as.numeric(x)

  stats::cor(log(x),log(y))^2
}

#' Compute the distance between two sets of vertices based on the MST
#' 
#' Specifically, focus on the largest connected component in the graph.
#' Return the number of edges so each vertex in \code{idx1} reaches its
#' closest \code{k} elements in \code{idx2}
#'
#' @param adj adjacency matrix
#' @param idx1 vector containing indices between 1 and \code{nrow(adj)}
#' @param idx2 vector containing indices between 1 and \code{nrow(adj)}
#' @param k positive integer
#'
#' @return numeric
#' @export
compute_mst_distance <- function(adj, idx1, idx2, k){
  g <- igraph::graph_from_adjacency_matrix(adj)
  tmp <- rep(1, igraph::vcount(g))
  group_vec <- rep(NA, igraph::vcount(g))
  group_vec[idx1] <- 1; group_vec[idx2] <- 2
  igraph::V(g)$group <- group_vec
  
  # find largest connected component
  tmp <- igraph::components(g)
  g <- igraph::induced_subgraph(g, v = which(tmp$membership == 1))
  g <- igraph::as.undirected(g)
  g <- igraph::minimum.spanning.tree(g)
  
  source <- which(igraph::V(g)$group == 1)
  sink <- which(igraph::V(g)$group == 2)
  dist_mat <- igraph::distances(g, v = source, to = sink)
  
  mean(apply(dist_mat, 1, function(x){sort(x, decreasing = F)[k]}))
}

#' Compute the distance between two sets of vertices based on graph root embedding
#'
#' @param eigenvectors eigenvectors of adjacency matrix
#' @param idx1 vector containing indices between 1 and \code{nrow(adj)}
#' @param idx2 vector containing indices between 1 and \code{nrow(adj)}
#' @param k number of positive and negative dimensions
#'
#' @return numeric
#' @export
compute_graph_root_distance <- function(eigenvectors, idx1, idx2, k){
  embedding_mat <- eigenvectors[,c(1:k, (ncol(eigenvectors)-k+1):ncol(eigenvectors))]
  
  dist_mat <- sapply(idx1, function(x){
    sapply(idx2, function(y){
      .l2norm(embedding_mat[x,] - embedding_mat[y,])
    })
  })
  
  mean(apply(dist_mat, 1, min))
}