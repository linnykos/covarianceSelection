rm(list=ls())
library(igraph)
load("../results/step5_res_alpha0.1.RData")
adj_our = adj_gene
g_our = igraph::graph_from_adjacency_matrix(adj_our, mode = "undirected")
aut_our = autism_genes

load("../results/step5_res_pfc35.RData")
adj_pfc = adj_gene
g_pfc = igraph::graph_from_adjacency_matrix(adj_pfc, mode = "undirected")
aut_pfc = autism_genes

load("../results/step5_res_alternative.RData")
adj_alt = adj_gene
g_alt = igraph::graph_from_adjacency_matrix(adj_alt, mode = "undirected")
aut_alt = autism_genes

#### let's plot some graphs?

# find all the autism genes with 2+ dnLoF or 1 dnLoF and is on iossifov

subset_graph <- function(g, aut){
  tada <- longitudinalGM::tada
  vec <- tada$Gene
  vec <- longitudinalGM::symbol_synonyms(vec)
  idx <- which(is.na(vec))
  tada <- tada[-idx,]; vec <- vec[-idx]
  tada$Gene <- vec

  genes_2dnlof <- intersect(tada$Gene[which(tada$dn.LoF >= 2)], aut)

  iossifov <- longitudinalGM::iossifov
  iossifov <- longitudinalGM::symbol_synonyms(iossifov)
  iossifov <- iossifov[!is.na(iossifov)]
  genes_1dnlof <- intersect(tada$Gene[which(tada$dn.LoF == 1)],
                            intersect(iossifov, aut))

  subset_genes <- sort(unique(c(genes_2dnlof, genes_1dnlof)))

  idx <- which(igraph::V(g)$name %in% subset_genes)
  dist_mat <- igraph::distances(g)
  all_idx <- sort(unique(unlist(lapply(idx, function(x){
    which(dist_mat[x,] <= 2)
  }))))

  gsub = igraph::induced_subgraph(g, all_idx)
  list(gsub = gsub, subset_genes = subset_genes)
}

gsub_our = subset_graph(g_our, aut_our)
gsub_pfc = subset_graph(g_pfc, aut_pfc)
gsub_alt = subset_graph(g_alt, aut_alt)

plot_graph <- function(g, aut, subset_genes){
  node_col <- rep("black", igraph::vcount(g))
  node_size <- rep(1, igraph::vcount(g))
  node_col[which(igraph::V(g)$name %in% aut)] <- rgb(205,40,54, max = 255)
  node_size[which(igraph::V(g)$name %in% aut)] <- 3
  node_size[which(igraph::V(g)$name %in% subset_genes)] <- 8

  set.seed(10)
  igraph::plot.igraph(g, vertex.size = node_size, vertex.label = NA,
                      vertex.color = node_col)
}

plot_graph(gsub_our$gsub, aut_our, gsub_our$subset_genes)
plot_graph(gsub_alt$gsub, aut_alt, gsub_alt$subset_genes)
plot_graph(gsub_pfc$gsub, aut_pfc, gsub_pfc$subset_genes)
