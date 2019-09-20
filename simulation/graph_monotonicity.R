
rm(list=ls())
###########################

# greedy triangle
greedy_triangle <- function(g){
  n <- igraph::vcount(g)
  g2 <- g
  
  h_seq <- vector("list", n-1)
  dens_vec <- rep(NA, length(h_seq))
  
  for(i in 1:length(h_seq)){
    h_seq[[i]] <- g2
    
    tri_mat <- matrix(igraph::triangles(g2), nrow=3)
    n_sub <- igraph::vcount(g2)
    dens_vec[[i]] <- ncol(tri_mat)/n_sub
      
    n_tri_count <- sapply(1:n_sub, function(i){length(which(tri_mat == i))})
    idx <- which.min(n_tri_count)
    
    g2 <- igraph::delete_vertices(g2, idx)
  }
  
  igraph::vcount(h_seq[[which.max(dens_vec)]])
}

greedy_triangle <- function(g){
  n <- igraph::vcount(g)
  g2 <- g
  
  h_seq <- vector("list", n-1)
  dens_vec <- rep(NA, length(h_seq))
  
  for(i in 1:length(h_seq)){
    h_seq[[i]] <- g2
    
    tri_mat <- matrix(igraph::triangles(g2), nrow=3)
    n_sub <- igraph::vcount(g2)
    dens_vec[[i]] <- ncol(tri_mat)/n_sub
    
    n_tri_count <- sapply(1:n_sub, function(i){length(which(tri_mat == i))})
    idx <- which.min(n_tri_count)
    
    g2 <- igraph::delete_vertices(g2, idx)
  }
  
  igraph::vcount(h_seq[[which.max(dens_vec)]])
}

random_deletion_edge <- function(g){
  while(TRUE){
    sample_idx <- sample(1:igraph::vcount(g), 1)
    neighbor_vec <- igraph::neighbors(g, sample_idx)
    
    if(length(neighbor_vec) > 0) break()
  }
  
  sample_idx2 <- as.numeric(sample(neighbor_vec, 1))
  
  igraph::delete_edges(g, igraph::get.edge.ids(g, c(sample_idx, sample_idx2)))
}

################################

set.seed(10)
n <- 50
combn_mat <- utils::combn(n, 2)
combn_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
g <- igraph::graph.empty(n = n, directed = F)
g <- igraph::add_edges(g, edges = combn_mat)
igraph::ecount(g)

trials <- 500
len_vec <- rep(NA, trials)
for(i in 1:trials){
  len_vec[i] <- greedy_triangle(g)
  g <- random_deletion_edge(g)
  print(paste0("Ecount: ", igraph::ecount(g), " for ", round( len_vec[i],2)))
}

set.seed(10)
n <- 50
combn_mat <- utils::combn(n, 2)
combn_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
g <- igraph::graph.empty(n = n, directed = F)
g <- igraph::add_edges(g, edges = combn_mat)
igraph::ecount(g)

trials <- 500
len_vec <- rep(NA, trials)
for(i in 1:trials){
  #len_vec[i] <- length(clique_selection(g)[[1]])
  g <- random_deletion_edge(g)
  print(paste0("Ecount: ", igraph::ecount(g), " for ", round( len_vec[i],2)))
}