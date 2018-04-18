rm(list=ls())
library(igraph)
load("../results/step5_res.RData")

# find all the autism genes with 2+ dnLoF or 1 dnLoF and is on iossifov
tada <- longitudinalGM::tada
vec <- tada$Gene
vec <- longitudinalGM::symbol_synonyms(vec)
idx <- which(is.na(vec))
tada <- tada[-idx,]; vec <- vec[-idx]
tada$Gene <- vec

genes_2dnlof <- intersect(tada$Gene[which(tada$dn.LoF >= 2)], autism_genes)

iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]
genes_1dnlof <- intersect(tada$Gene[which(tada$dn.LoF == 1)],
                          intersect(iossifov, autism_genes))

subset_genes <- sort(unique(c(genes_2dnlof, genes_1dnlof)))

# find the graph
idx <- which(colnames(adj_gene) %in% subset_genes)
g <- igraph::graph_from_adjacency_matrix(adj_gene, mode = "undirected")
dist_mat <- igraph::distances(g)
all_idx <- sort(unique(unlist(lapply(idx, function(x){
  which(dist_mat[x,] <= 2)
}))))

#take subgraph
g2 <- igraph::induced_subgraph(g, all_idx)

#find largest connected component and take subgraph again
comp <- igraph::components(g2)
nam <- names(which(comp$membership == which.max(comp$csize)))
g3 <- igraph::induced_subgraph(g2, nam)

igraph::vcount(g3)/igraph::ecount(g3)
igraph::vcount(g3)/(igraph::ecount(g3)*(igraph::ecount(g3)-1)/2)
length(which(igraph::V(g3)$name %in% autism_genes))/igraph::ecount(g3)

node_col <- rep("black", igraph::vcount(g3))
node_size <- rep(1, igraph::vcount(g3))
node_col[which(igraph::V(g3)$name %in% autism_genes)] <- rgb(205,40,54, max = 255)
node_size[which(igraph::V(g3)$name %in% autism_genes)] <- 3
node_size[which(igraph::V(g3)$name %in% subset_genes)] <- 8

pdf("../figures/figure_12.pdf", height = 5, width = 5)
set.seed(10)
igraph::plot.igraph(g3, vertex.size = node_size, vertex.label = NA,
                    vertex.color = node_col)
graphics.off()

##################

# compute the distances among
idx <- which(colnames(adj_gene) %in% subset_genes)
g <- igraph::graph_from_adjacency_matrix(adj_gene, mode = "undirected")
dist_mat <- igraph::distances(g)
dist_mat <- dist_mat[idx, idx]
row_keep <- which(apply(dist_mat, 1, function(x){sum(is.infinite(x)) < length(x)-1}))
dist_mat <- dist_mat[row_keep, row_keep]
median(dist_mat[lower.tri(dist_mat, diag = F)])
