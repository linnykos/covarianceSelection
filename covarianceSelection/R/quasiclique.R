tsourakakis_2013 <- function(g, threshold = 0.95, iter_max = round(igraph::vcount(g)/2)){
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
  
  as.numeric(node_set)
}

# if node_candidate = NA, we are subtracting, not adding
.find_candidate <- function(g, threshold, node_set, node_candidate, den_org){
  stopifnot(length(node_set) > 1, is.character(node_set))
  
  if(all(is.na(node_candidate))){
    node_set <- sample(node_set)
    n2 <- length(node_set)
  } else {
    stopifnot(!any(node_set %in% node_candidate))
    node_candidate <- sample(node_candidate)
    n2 <- length(node_candidate)
  }
  
  for(i in 1:n2){
    if(all(is.na(node_candidate))){
      next_set <- node_set[-i]
    } else {
      next_set <- c(node_set, node_candidate[i])
    }
    
    den_new <- .tsourakakis_obj(g, threshold, next_set)
    if(den_new > den_org) break()
  }
  
  if(den_new <= den_org) return(NA)
  sort(next_set)
}

.tsourakakis_obj <- function(g, threshold, node_set){
  stopifnot(is.character(node_set))
  
  g_sub <- igraph::induced_subgraph(g, node_set)
  igraph::ecount(g_sub) - threshold * choose(length(node_set), 2)
}

.tsourakakis_initialize <- function(g){
  tri_mat <- matrix(igraph::triangles(g), nrow=3)
  n <- igraph::vcount(g)
  n_tri_count <- sapply(1:n, function(i){length(which(tri_mat == i))})
  deg_vec <- igraph::degree(g)
  
  idx <- which.max(n_tri_count/deg_vec)
  as.character(igraph::V(g)$name[idx])
}

########################

chen_2010 <- function(g, threshold = 0.95){
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
  
  as.numeric(node_set_internal)
}

# node set should be characters for safety
.chen_check_density <- function(g, node_set){
  stopifnot(is.character(node_set), all(node_set %in% as.character(igraph::V(g)$name)))
  if(length(node_set) == 1) return(0)
  
  g_sub <- igraph::induced_subgraph(g, node_set)
  igraph::ecount(g_sub)/choose(igraph::vcount(g_sub), 2)
}

.chen_separate <- function(g, c_matrix, threshold = 0.95, check = F){
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
  
  list(g1 = g1, g2 = g2, c_matrix1 = c_matrix1, c_matrix2 = c_matrix2)
}

.chen_object <- function(g, c_matrix, check = F){
  node_set <- as.numeric(igraph::V(g)$name)
  
  if(check) stopifnot(all(c_matrix[,1] %in% node_set), all(c_matrix[,2] %in% node_set))
  
  structure(list(g = g, node_set = node_set, c_matrix = c_matrix),
            class = "chen_object")
}

# diagonal of adj must be 1 to prevent the computed quantity from being undefined
.form_c_matrix <- function(g){
  adj <- as.matrix(igraph::as_adj(g))
  diag(adj) <- 1
  n <- nrow(adj)
  
  combn_mat <- utils::combn(n, 2)
  c_matrix <- matrix(0, nrow = ncol(combn_mat), ncol = 3)
  c_matrix[,1:2] <- t(combn_mat)
  for(i in 1:nrow(c_matrix)){
    idx1 <- c_matrix[i,1]; idx2 <- c_matrix[i,2]
    c_matrix[i,3] <- adj[idx1,]%*%adj[idx2,]/(.l2norm(adj[idx1,])*.l2norm(adj[idx2,]))
  }
  
  c_matrix[order(c_matrix[,3]),,drop = F]
}

########

anderson_2009 <- function(g){
  g <- igraph::as.undirected(g)
  g <- igraph::simplify(g)
  n <- igraph::vcount(g)
  igraph::V(g)$name <- 1:n
  
  h_seq <- vector("list", n-1)
  dens_vec <- rep(NA, length(h_seq))
  h_seq[[1]] <- g
  
  for(i in 1:length(h_seq)){
    total_deg <- igraph::ecount(h_seq[[i]])
    n_sub <- igraph::vcount(h_seq[[i]])
    dens_vec[[i]] <- total_deg/n_sub
    
    deg_vec <- igraph::degree(h_seq[[i]])
    idx <- which.min(deg_vec)
    
    h_seq[[i+1]] <- igraph::delete_vertices(h_seq[[i]], idx)
  }
  
  as.numeric(igraph::V(h_seq[[which.max(dens_vec)]])$name)
}

#####

tsourakakis_2014_exact <- function(g){
  g <- igraph::as.undirected(g)
  g <- igraph::simplify(g)
  
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

######

tsourakakis_2014_approximate <- function(g){
  g <- igraph::as.undirected(g)
  g <- igraph::simplify(g)
  n <- igraph::vcount(g)
  igraph::V(g)$name <- 1:n
  
  h_seq <- vector("list", n-2)
  dens_vec <- rep(NA, length(h_seq))
  h_seq[[1]] <- g
  
  for(i in 1:length(h_seq)){
    tri_mat <- matrix(igraph::triangles(h_seq[[i]]), nrow=3)
    n_sub <- igraph::vcount(h_seq[[i]])
    dens_vec[[i]] <- ncol(tri_mat)/n_sub
    
    n_tri_count <- sapply(1:n_sub, function(i){length(which(tri_mat == i))})
    idx <- which.min(n_tri_count)
    
    h_seq[[i+1]] <- igraph::delete_vertices(h_seq[[i]], idx)
  }
  
  as.numeric(igraph::V(h_seq[[which.max(dens_vec)]])$name)
}