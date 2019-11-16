if(verbose) print(paste0(Sys.time(), "Start of step 8: Compiling the results"))

validated_genes <- covarianceSelection::validated_genes$Gene

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) 
num_our <- length(intersect(genes_our, validated_genes)) 
num_intersect <- length(intersect(genes_our_intersect, validated_genes)) 

# output some summaries
c(num_pfc35, num_nodawn, num_our, num_intersect)
c(length(genes_pfc35), length(genes_nodawn), length(genes_our), length(genes_our_intersect))
c(sum(as.matrix(adj_pfc35))/2, sum(as.matrix(adj_our))/2)
c(covarianceSelection::compute_scale_free(as.matrix(adj_pfc35)), covarianceSelection::compute_scale_free(as.matrix(adj_our)))

cbind(names(dat_list[idx_our]), sapply(dat_list[idx_our], nrow))
covarianceSelection::binning(names(dat_list)[idx_our])

# computing eigen decompositions for the figure diagnostics
eigen_pfc35 <- eigen(as.matrix(adj_pfc35))
eigen_our <- eigen(as.matrix(adj_our))

save.image(file = paste0(save_filepath, "/step8_results", filepath_suffix, ".RData"))

