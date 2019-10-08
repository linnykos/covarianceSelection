if(verbose) print(paste0(Sys.time(), "Start of step 7: Compiling the results"))
validated_genes <- covarianceSelection::validated_genes$Gene

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  # 6 genes it seems?
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) # it seems to be 32 genes?
num_all <- length(intersect(genes_all, validated_genes)) # it seems to be 8 genes?
num_our <- length(intersect(genes_our, validated_genes)) # seems to be 17 genes?

########

