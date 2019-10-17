if(verbose) print(paste0(Sys.time(), "Start of step 7: Compiling the results"))
validated_genes <- covarianceSelection::validated_genes$Gene

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  # 28/175 genes?
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) # 11/13 genes?
num_all <- length(intersect(genes_all, validated_genes)) # 28/102 genes?
num_our <- length(intersect(genes_our, validated_genes)) # 30/179 genes?

c(num_pfc35, num_nodawn, num_all, num_our)
c(length(genes_pfc35), length(genes_nodawn), length(genes_all), length(genes_our))
c(hmrf_pfc35$c, hmrf_all$c, hmrf_our$c)
c(sum(as.matrix(adj_pfc35))/2, sum(as.matrix(adj_all))/2, sum(as.matrix(adj_our))/2)

cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_pfc35, scale_vec_pfc35)
cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_all, scale_vec_all)
cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_our, scale_vec_our)

rm(list = c("validated_genes"))

save.image(file = paste0(save_filepath, "/step7_results", filepath_suffix, ".RData"))

