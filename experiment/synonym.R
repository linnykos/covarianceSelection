rm(list=ls())
load("../results/step5_res_alternative.RData")
zz <- as.character(autism_genes)
zreport <- report
zg <- igraph::graph_from_adjacency_matrix(adj_gene)
load("../results/step5_res_pfc35.RData")
yy <- as.character(autism_genes)
yreport <- report
yg <- igraph::graph_from_adjacency_matrix(adj_gene)
load("../results/step5_res.RData")
xx <- as.character(autism_genes)
xreport <- report
xg <- igraph::graph_from_adjacency_matrix(adj_gene)

iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

length(intersect(xx, iossifov))
length(intersect(yy, iossifov))
length(intersect(zz, iossifov))

length(intersect(xx, yy))
length(intersect(xx, zz))
length(intersect(yy, zz))

########################

# candidates to swap out
candidates_out <- xx[which(!xx %in% c(yy, zz))]
candidates_out <- candidates_out[!candidates_out %in% iossifov]

# create a string that can be pasted
# print(paste0(candidates_out, collapse = "\',\'"))
# idx <- which(colnames(genexp) %in% )
# cbind(colnames(genexp)[idx], idx)

# check to see the candidates are not too high on the report
report_list <- list(zreport, yreport, xreport)
for(i in 1:3){
  report_list[[i]] <- report_list[[i]][order(report_list[[i]]$FDR),]
  vec <- which(report_list[[i]]$Gene %in% candidates_out)
  names(vec) <- report_list[[i]]$Gene[vec]
  print(vec)
  print("------")
}
# see how close the canddiates are to the others in desired network
swap_out <- ?????
dist_mat <- igraph::distances(xg)
dist_mat <- dist_mat[which(rownames(dist_mat) %in% candidates_out), which(colnames(dist_mat) %in% intersect(iossifov, xx))]
thre_vec <- apply(dist_mat, 1, function(x){length(which(x <= 3))})

###########

# start to look for candidates to swap in:
candidates_in <- iossifov[!iossifov %in% c(zz, yy, xx)]

# take in candidates already found
liliu <- read.csv("data-raw/aoas844_supp.csv", sep = ";")
liliu <- liliu$Gene

candidates_in[!candidates_in %in% liliu]

# take candidates with a dnlof already
tada <- longitudinalGM::tada
candidates_in[which(candidates_in %in% tada[which(tada$dn.LoF >= 1),]$Gene)]

tada <- tada[which(tada$Gene %in% candidates_in),]
tada <- tada[order(tada$pval.TADA),]

# check positions


# check to see if swap_in are not in graph?
swap_in <- c("DIP2A", "MED13L", "RIMS1") # positions: 4853, 3551, 6085
genes_graph <- colnames(igraph::distances(xg))
which(swap_in %in% genes_graph)

# check to see swap_in is irrelevant in report
xreport <- xreport[order(xreport$FDR),]
xreport[which(xreport$Gene %in% swap_in),]
xreport[which(xreport$Gene %in% intersect(xx, iossifov)),]
# currently in xx and iossifov:
c("ANK2","DYRK1A","ARID1B","CHD8","ADNP","POGZ","WDFY3","NCKAP1","RANBP17",
  "KDM5B","KDM6B","FOXP1","DSCAM","SPAST","MED13L","PHF2")

# mapping:
# pretty good: got 18 vs 16 vs 16
# DICER1 -> MED13L
# ITGAV -> DIP2A
# ZMYND11 -> RIMS1

##########
# pretty good: got 19 vs 16 (pfc) vs 17
REPIN1 -> MED13L
ITGAV -> DIP2A
ZMYND11 -> RIMS1
