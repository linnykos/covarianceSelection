#' Spectral clustering
#'
#' @param g \code{igraph} object
#' @param K number of clusters
#'
#' @return vector
#'
#' @reference
#' \url{https://github.com/lingxuez/splitSpectral/blob/master/splitSpectral.R}
#' Lei, Jing, and Lingxue Zhu. "Generic Sample Splitting for Refined Community Recovery in Degree Corrected Stochastic Block Models." (2016).
.spectral_cluster <- function(g, K = 2){
  adj <- igraph::as_adj(g)
  diag(adj) <- 0

  U <- eigen(adj)$vectors[,1:K]
  U <- .row_norm(U)
  stats::kmeans(U, centers=K, nstart=20)$cluster
}


#' Row normalization
#'
#' Normalize the non-zero rows of a matrix A
#' to have unit L2 norm
#'
#' @param A the input matrix
#'
#' @return the new matrix after row-normalization.
.row_norm <- function(A) {
  A_rownorm <- apply(A, 1, .l2norm)
  i_nonzero <- which(A_rownorm > 0)
  A[i_nonzero, ] <- apply(A[i_nonzero, ], 2, function(x){x/.l2norm(x)})
  A
}

#' Find largest clustering based on spectral clustering
#'
#' Wrapper for \code{.spectral_cluster}
#'
#' @param g \code{igraph} object
#' @param K_vec vector of number of clusters to consider
#' @param threshold number between 0 and 1 for density of quasi-clique
#'
#' @return vector
#' @export
spectral_selection <- function(g, K_vec = 2:5, threshold = 0.95){
  res <- sapply(K_vec, function(K){
    clustering <- .spectral_cluster(g, K)
    
    size_vec <- table(clustering)
    den_vec <- sapply(1:K, function(i){
      idx <- which(clustering == i)
      igraph::ecount(igraph::induced_subgraph(g, idx))/(choose(length(idx), 2))
    })
    
    idx <- which(den_vec >= threshold)
    if(length(idx) == 0) return(numeric(0))
    cluster_idx <- (c(1:K)[idx])[which.max(size_vec[idx])]
    
    which(clustering == cluster_idx)
  })
  
  res[[which.max(sapply(res, length))]]
}
