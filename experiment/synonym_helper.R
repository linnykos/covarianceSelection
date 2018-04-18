rm(list=ls())
iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

load("step5_res.RData")
genes1 <- autism_genes
adj1 <- adj_gene

load("step5_res_pfc35.RData")
genes2 <- autism_genes
adj2 <- adj_gene

tada <- longitudinalGM::tada
vec <- tada$Gene
vec <- longitudinalGM::symbol_synonyms(vec)
idx <- which(is.na(vec))
tada <- tada[-idx,]; vec <- vec[-idx]
tada$Gene <- vec

#possible genes
working_genes <- sort(setdiff(genes1, genes2))

#swap genes in
genes1plus <- intersect(tada$Gene[which(tada$dn.LoF >= 1)], iossifov)
genestarget <- setdiff(genes1plus, c(as.character(genes1), as.character(genes2)))
#DIP2A, ILF2, RIMS1

#include gene already in pfc35
genestarget <- c(genestarget, intersect(setdiff(genes2, genes1), iossifov))
#NINL

#include 2 most promisin genes in TADA not already included
iossifov_remaining <- setdiff(iossifov, genestarget)
tada_pval <- unlist(sapply(iossifov_remaining,
                              function(x){tada[which(tada$Gene == x),"pval.TADA"]}))
genestarget <- c(genestarget, iossifov_remaining[which(rank(tada_pval, ties.method = "random") <= 2)])
#ADNP, ANK2

load("/raid6/Kevin/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)

#determine brain-expressed genes
brain_expression <- longitudinalGM::brain_expression
brain_genes <- brain_expression$Gene[brain_expression$Brain_expressed != 'No']
idx <- which(colnames(genexp) %in% brain_genes)
genexp <- genexp[,idx]
gene_order <- colnames(genexp)

rbind(gene_order[which(gene_order %in% genestarget)], which(gene_order %in% genestarget))
rbind(gene_order[which(gene_order %in% working_genes)], which(gene_order %in% working_genes))

#find genes with no neighbors
working_genes2 <- gene_order[which(gene_order %in% working_genes)][1:25]

num_neighbors <- matrix(0, ncol = length(working_genes2), nrow = 2)
for(i in 1:length(working_genes2)){
  num_neighbors[1,i] <- sum(adj1[which(rownames(adj1) == working_genes2[i]),])
  num_neighbors[2,i] <- sum(adj2[which(rownames(adj2) == working_genes2[i]),])
}
colnames(num_neighbors) <- working_genes2
#NCOR1, NCOA3, ZMYND11, DEAF1, SLCO6A1, DYNC1H1

#mapping:
#NCOR1 -> ANK2
#DYNC1H1 -> DIP2A
#DEAF1 -> ILF2
#NCOA3 -> NINL
#ZMYND11 -> RIMS1
#SLCO6A1 -> ADNP
