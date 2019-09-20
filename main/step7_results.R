validated_genes <- read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))
validated_genes <- covarianceSelection::symbol_synonyms(validated_genes, verbose = T)

length(intersect(genes_pfc35, validated_genes)) #16 genes it seems?
# without dawn, it seems to be 30 genes?