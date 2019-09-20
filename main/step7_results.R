validated_genes <- read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))
validated_genes <- covarianceSelection::symbol_synonyms(validated_genes, verbose = T)

length(intersect(genes_pfc35, validated_genes))  # 5 genes it seems?
length(intersect(genes_nodawn, validated_genes)) # it seems to be 32 genes?
