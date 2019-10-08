## code to prepare `validated_genes` dataset goes here
rm(list=ls())
validated_genes <- read.csv("data-raw/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))
validated_genes <- covarianceSelection::symbol_synonyms(validated_genes, verbose = T)
validated_genes <- as.data.frame(validated_genes)
colnames(validated_genes) <- "Gene"

usethis::use_data(validated_genes, overwrite = T)
