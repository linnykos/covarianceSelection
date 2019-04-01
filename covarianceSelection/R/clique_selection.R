#' Find largest soft clique in graph
#'
#' For a given graph \code{g} (igraph object) and a threshold \code{threshold}
#' (number between 0 and 1), find the largest set of nodes such that
#' the number of edges among this set of nodes in \code{g} is at least
#' \code{threshold} percentage of a full clique.
#'
#' This is returned as a list of indices, so if there are multiple largest sets
#' of the same size, the list will have multiple set of indices.
#' 
#' @param g \code{igraph} object
#' @param num_pos numeric vector
#' @param num_neg numeric vector
#' @param threshold numeric
#'
#' @return list of indices
#' @export
clique_selection <- function(g, num_pos = c(1:4), num_neg = c(0:3), threshold = 0.95){
  combn_mat <- as.matrix(expand.grid(num_pos, num_neg))
  stopifnot(apply(combn_mat, 1, function(x){any(x > 0)}))
  
  res <- lapply(1:nrow(combn_mat), function(x){
    tryCatch({
      covarianceSelection:::.clique_selection(g, combn_mat[x,1], combn_mat[x,2], threshold)
    }, error = function(e){
      list(numeric(0))
    })
  })
  
  #unravel the list
  k <- 1
  res2 <- vector("list", 0)
  for(i in 1:length(res)){
    if(length(res[[i]]) > 0){
      for(j in 1:length(res[[i]])){
        res2[[k]] <- res[[i]][[j]]
        k <- k+1
      }
    }
  }
  
  len_vec <- sapply(res2, length)
  res2[which(len_vec == max(len_vec))]
}

#' Select clique based on intersection
#'
#' @param lis list of indices
#' @param idx vector of indices
#' @param adj adjacency matrix
#'
#' @return one element of \code{lis}
#' @export
select_clique <- function(lis, idx, adj){
  stopifnot(max(idx) <= ncol(adj), ncol(adj) == nrow(adj))
  stopifnot(all(sapply(lis, max) <= ncol(adj)))
  
  if(length(lis) > 1){
    vec <- sapply(lis, function(x){
      length(which(idx %in% x))
    })
    
    if(length(which(vec == max(vec))) > 1){
      idx <- which.max(sapply(lis, function(x){sum(adj[x])}))
      lis[[idx]]
    } else {
      lis[[which(vec == max(vec))]]
    }
  } else {
    lis[[1]]
  }
}

######################

.clique_selection <- function(g, num_pos = 2, num_neg = 0, threshold = 0.95){
  adj_mat <- as.matrix(igraph::as_adjacency_matrix(g))
  diag(adj_mat) <- 0
  
  eigen_res <- eigen(adj_mat)
  pos_mat <- .form_eigen_mat(eigen_res, num = num_pos, pos = TRUE)
  neg_mat <- .form_eigen_mat(eigen_res, num = num_neg, pos = FALSE)
  
  dist_mat <- stats::as.dist(.form_dist_mat(pos_mat, neg_mat))
  hclust_res <- stats::hclust(dist_mat)
  lis <- .extract_level_set(hclust_res)
  stopifnot(length(lis) >= 1)
  
  bool_vec <- sapply(lis, function(x){
    .pass_threshold(adj_mat[x,x], threshold)
  })
  
  size_vec <- sapply(lis, length)
  names(size_vec) <- 1:length(lis)
  size_vec <- size_vec[bool_vec]
  
  max_val <- max(size_vec)
  idx <- which(size_vec == max_val)
  
  res <- lapply(idx, function(x){
    sort(lis[[as.numeric(names(size_vec[x]))]])
  })
  names(res) <- NULL
  
  res
}

.form_eigen_mat <- function(eigen_res, num, pos = TRUE){
  n <- nrow(eigen_res$vectors)
  if(num == 0){
    return(matrix(0, nrow = n, ncol = 1))
  }
  
  if(pos){
    idx <-  which(eigen_res$values > 0)
  } else {
    idx <- rev(which(eigen_res$values < 0))
  }
  
  stopifnot(length(idx) >= num)
  
  if(num == 1){
    diag_mat <- matrix(sqrt(abs(eigen_res$values[idx[1]])), 1, 1)
  } else {
    diag_mat <- diag(sqrt(abs(eigen_res$values[idx[1:num]])))
  }
  
  eigen_res$vectors[,idx[1:num],drop = F] %*% diag_mat
}

.form_dist_mat <- function(pos_mat, neg_mat){
  stopifnot(nrow(pos_mat) == nrow(neg_mat))
  
  n <- nrow(pos_mat)
  
  mat <- matrix(0, n, n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      mat[i,j] <- .krein_dist_func(pos_mat[i,], neg_mat[i,], pos_mat[j,], neg_mat[j,])
      mat[j,i] <- mat[i,j]
    }
  }
  
  mat
}

.krein_dist_func <- function(pos_vec1, neg_vec1, pos_vec2, neg_vec2){
  stopifnot(length(pos_vec1) == length(pos_vec2),
            length(neg_vec1) == length(neg_vec2))
  
  sqrt(.l2norm(pos_vec1 - pos_vec2)^2 + .l2norm(neg_vec1 - neg_vec2)^2)
}

.extract_level_set <- function(hclust_res){
  lis <- vector("list", 0)
  len <- 1
  
  m <- nrow(hclust_res$merge)
  for(i in 1:m){
    idx <- hclust_res$merge[i,]
    if(all(idx < 0)){
      lis[[len]] <- abs(idx)
    } else {
      if(idx[1] < 0) vec1 <- abs(idx[1]) else vec1 <- lis[[idx[1]]]
      if(idx[2] < 0) vec2 <- abs(idx[2]) else vec2 <- lis[[idx[2]]]
      lis[[len]] <- c(vec1, vec2)
    }
    
    len <- len + 1
  }
  
  lis
}

.pass_threshold <- function(adj_mat, threshold){
  n <- nrow(adj_mat)
  num_threshold <- threshold*n*(n-1)/2
  
  sum(adj_mat)/2 > num_threshold
}
