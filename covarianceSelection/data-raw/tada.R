## code to prepare `validated_genes` dataset goes here
rm(list=ls())
tada <- read.csv("data-raw/TADA_Results_2231trios_1333trans_1601cases_5397controls_March26_pvalues.csv") # 18735 genes
tada <- tada[,which(colnames(tada) %in% c("Gene", "dn.LoF", "qvalue", "pval.TADA"))]
vec <- covarianceSelection::symbol_synonyms( tada$Gene, verbose = T)
unknown_genes_idx <- which(sapply(vec, length) == 0)
tada <- tada[-unknown_genes_idx,] # 18700 genes
vec <- vec[-unknown_genes_idx]; vec <- unlist(vec)
tada$Gene <- vec

#remove duplicated tada by keeping the one with the lowest p-value
tada <- tada[-which(duplicated(tada$Gene)),] #18498 genes
usethis::use_data(tada, overwrite = T)
