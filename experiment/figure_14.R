rm(list=ls())
dist_list <- vector("list", 3)
additional_name_vec <- c("", "_alternative", "_pfc35")

for(counter in 1:3){
  load(paste0("../results/step5_res", additional_name_vec[counter], ".RData"))

  idx <- which(colnames(adj_gene) %in% autism_genes)
  g <- igraph::graph_from_adjacency_matrix(adj_gene, mode = "undirected")
  dist_mat <- igraph::distances(g)
  dist_mat <- dist_mat[idx, idx]
  raw_vec <- table(dist_mat[upper.tri(dist_mat, diag = F)])

  #now simulate the base
  set.seed(1)
  g_sim <- igraph::sample_gnm(igraph::vcount(g), igraph::ecount(g))
  dist_mat <- igraph::distances(g)
  dist_mat <- dist_mat[1:length(autism_genes), 1:length(autism_genes)]
  normalization <- table(dist_mat[upper.tri(dist_mat, diag = F)])

  val <- as.numeric(names(raw_vec))
  vec <- sort(unlist(lapply(1:length(val), function(x){
    denom <- sum(normalization)
    count <- raw_vec[x]
    lower_idx <- which(as.numeric(names(normalization)) < val[x])
    base_num <- ifelse(length(lower_idx) == 0, 0, sum(normalization[lower_idx]))

    exact_idx <- which(as.numeric(names(normalization)) == val[x])
    if(length(exact_idx) == 1){
      numer <- base_num + seq(1, normalization[exact_idx], length.out = count)
    } else {
      numer <- rep(base_num, length = count)
    }

    numer/denom
  })))

  dist_list[[counter]] <- vec
}

hist(dist_list[[1]], col = "gray")
hist(dist_list[[2]], col = rgb(1, 0, 0, 0.5), add = T)
hist(dist_list[[3]], col = rgb(0, 1, 0, 0.5), add = T)
