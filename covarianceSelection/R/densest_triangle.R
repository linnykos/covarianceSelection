densest_triangle <- function(g){
  stopifnot(igraph::vcount(g) > 0, igraph::ecount(g) > 0)
  n <- igraph::vcount(g)
  tri_mat <- matrix(igraph::triangles(g), nrow=3)
  n_tri <- ncol(tri_mat)
  s_idx <- n+n_tri+1
  t_idx <- s_idx+1
  
  l <- 0; u <- n^3
  iter <- 1
  
  while(u >= l + 1/(n*(n-1))){
    print(iter)
    alpha <- (l+u)/2
    n_tri_count <- sapply(1:n, function(i){length(which(tri_mat == i))})
    
    # construct the adjacency matrix
    adj_mat <- matrix(0, t_idx, t_idx)
    for(i in 1:ncol(tri_mat)){
      for(j in 1:nrow(tri_mat)){
        idx1 <- tri_mat[j,i] # idx of node
        idx2 <- n+i # idx of triangle
        adj_mat[idx1, idx2] <- 1
        adj_mat[idx2, idx1] <- 2
      }
    }
    
    for(i in 1:n){
      adj_mat[s_idx, i] <- n_tri_count[i]
    }
    
    for(i in 1:n_tri){
      adj_mat[i+n, t_idx] <- 3*alpha
    }
    
    ig <- igraph::graph.adjacency(adj_mat, mode="directed", weighted=TRUE)
    res <- igraph::st_min_cuts(ig, source = s_idx, target = t_idx)
    
    # keep the smallest partition
    idx <- which.min(sapply(res$partition1s, length))
    partition <- as.numeric(res$partition1s[[idx]])
    #print(partition)
    
    if(all(partition == s_idx)) {
      u <- alpha
    } else {
      l <- alpha
      set <- intersect(partition[-which(partition == s_idx)], 1:n)
      print(set)
    }
    
    iter <- iter + 1
  }
  
  set
}