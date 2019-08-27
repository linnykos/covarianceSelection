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

  U <- mgcv::slanczos(adj, K)$vectors
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
  A_rownorm <- sqrt(rowSums(A^2))
  i_nonzero <- which(A_rownorm > 0)
  A[i_nonzero, ] <- apply(A[i_nonzero, ], 2, function(x){x / A_rownorm[i_nonzero]})
  A
}

#' Find largest clustering based on spectral clustering
#'
#' Wrapper for \code{.spectral_cluster}
#'
#' @param g \code{igraph} object
#' @param K number of clusters
#' @param target_idx which indices to select cluster based on
#'
#' @return vector
#' @export
spectral_selection <- function(g, K = 2, target_idx = NA){
  clustering <- .spectral_cluster(g, K)
  if(any(!is.na(target_idx))){
    tab <- table(clustering[target_idx])
  } else {
    tab <- table(clustering)
  }

  idx <- as.numeric(names(tab)[which.max(tab)])
  which(clustering == idx)
}
