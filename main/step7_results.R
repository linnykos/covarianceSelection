if(verbose) print(paste0(Sys.time(), "Start of step 7: Compiling the results"))

validated_genes <- read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))
validated_genes <- covarianceSelection::symbol_synonyms(validated_genes, verbose = T)

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  # 32 genes it seems?
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) # it seems to be 32 genes?
num_all <- length(intersect(genes_all, validated_genes)) # it seems to be 14 genes?
num_our <- length(intersect(genes_our, validated_genes)) # seems to be 25 genes?


