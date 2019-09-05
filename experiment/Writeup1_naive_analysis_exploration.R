rm(list=ls())
load("../results/step2_naive_analysis.RData")
validated_genes  = read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes = as.character(as.matrix(validated_genes))

length(intersect(autism_genes, validated_genes))
gene_names <- sort(intersect(autism_genes, validated_genes))
tada[which(tada$Gene %in% gene_names),]

tada_full <- read.csv("../../raw_data/TADA_Results_2231trios_1333trans_1601cases_5397controls_March26_pvalues.csv")
tada_full[which(tada_full$Gene %in% gene_names),]

zz = tada_full[which(tada_full$Gene %in% autism_genes),]
table(zz$dn.LoF)

neigh_vec <- sapply(autism_genes, function(x){
  idx <- which(tada$Gene == x)
  neigh_idx <- which(adj_gene[idx,] != 0)
  tada_idx <- which(tada_full$Gene %in% names(neigh_idx))
  length(which(tada_full[tada_idx, "dn.LoF"] >= 1))
})

p_val_vec <- sapply(autism_genes, function(x){
  tada[which(tada$Gene == x), "pval.TADA"]
})
plot(p_val_vec, neigh_vec)

p_val_after <- sapply(autism_genes, function(x){
  report[which(report$Gene == x), "FDR"]
})

plot(p_val_vec, p_val_after)
