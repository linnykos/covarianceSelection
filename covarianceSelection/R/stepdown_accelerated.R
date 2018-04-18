.compute_max_accelerated <- function(num_list, combn_mat){
  g <- .construct_graph(combn_mat, value = Inf)

  #construct spanning tree (shuffle the indices and make a chain)
  tree <- igraph::mst(g)
  tree_vertex_mat <- .return_vertices(tree)

  #create the list of edges to be formed
  diff_g <- igraph::difference(g, tree)
  vertex_mat <- .return_vertices(diff_g)

  #compute initial set of bootstrap stat (so for this, we need all the ingredients to compute the numerator)
  val <- .compute_difference_from_edges(num_list, tree_vertex_mat)
  g <- .set_edge_weight(g, tree_vertex_mat, val)

  #record current max
  max_val <- max(abs(val))

  #while
  for(i in 1:ncol(vertex_mat)){
    #compute shortest path
    val <- igraph::distances(g, v = vertex_mat[1,i], to = vertex_mat[2,i])

    #compare to maximum value
    if(val > max_val){
      #compute its test statistic
      val <- .compute_covStat(num_list[[as.numeric(vertex_mat[1,i])]],
                              num_list[[as.numeric(vertex_mat[2,i])]], 1, 1, squared = F)
      g[vertex_mat[1,i], vertex_mat[2,i]] <- val
      g[vertex_mat[2,i], vertex_mat[1,i]] <- val

      #see if its larger than the current max
      if(val > max_val) max_val <- val
    }
  }

  #return
  as.numeric(max_val)
}

.construct_graph <- function(edge_mat, value = 1){
  if(length(edge_mat) == 1) edge_mat <- utils::combn(edge_mat, 2)

  idx <- as.character(sort(unique(unlist(as.vector(edge_mat)))))
  n <- length(idx)
  g <- igraph::graph.empty(n = n, directed = F)
  igraph::V(g)$name <- idx
  g <- igraph::add_edges(g, edges = t(apply(edge_mat, 1, as.character))) #is this better than full dense?
  igraph::set_edge_attr(g, "weight", value = value)
}

.return_vertices <- function(g){
  edgelist <- igraph::E(g)
  sapply(edgelist, function(x){
    as.character(range(as.numeric(igraph::ends(g, x))))
  })
}

.compute_difference_from_edges <- function(num_list, vertex_mat){
  vertex_mat <- apply(vertex_mat, 2, as.numeric)

  apply(vertex_mat, 2, function(x){
    .compute_covStat(num_list[[x[1]]], num_list[[x[2]]], 1, 1, squared = F)
  })
}

.set_edge_weight <- function(g, vertex_mat, val){
  stopifnot(ncol(vertex_mat) == length(val))
  for(i in 1:ncol(vertex_mat)){
    g[vertex_mat[1,i], vertex_mat[2,i]] <- val[i]
    g[vertex_mat[2,i], vertex_mat[1,i]] <- val[i]
  }

  g
}
