.remove_edges_uniform <- function(g, remove_num){
  idx <- sample(1:igraph::ecount(g), remove_num)
  igraph::delete_edges(g, idx)
}

.remove_edges_proportional <- function(g, remove_num){
  n <- igraph::vcount(g)
  for(i in 1:remove_num){
    deg_vec <- igraph::degree(g)
    
    prob_vec <- 1/deg_vec
    prob_vec[is.infinite(prob_vec)] <- 0
    
    idx <- sample(1:n, 1, prob = prob_vec)
    
    # find neighbors
    neigh_idx <- igraph::neighbors(g, idx)
    if(length(neigh_idx) == 1){
      neigh_selected = neigh_idx
    } else{
      neigh_selected <- sample(neigh_idx, 1)
    }
    
    g <- igraph::delete_edges(g, paste0(idx, "|", neigh_selected))
  }
  
  g
}
