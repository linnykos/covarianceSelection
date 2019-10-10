if(verbose) print(paste0(Sys.time(), "Start of step 7: Compiling the results"))
validated_genes <- covarianceSelection::validated_genes$Gene

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  # 32 genes it seems?
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) # it seems to be 32 genes?
num_all <- length(intersect(genes_all, validated_genes)) # it seems to be 34 genes?
num_our <- length(intersect(genes_our, validated_genes)) # seems to be 31 genes?

########

# tmp_pfc35 <- as.character(report_pfc35[which(report_pfc35$FDR <= 0.01), "Gene"])
# tmp_all <- as.character(report_all[which(report_all$FDR <= 0.01), "Gene"])
# tmp_our <- as.character(report_our[which(report_our$FDR <= 0.01), "Gene"])
# 
# length(intersect(tmp_pfc35, validated_genes))
# length(intersect(tmp_all, validated_genes))
# length(intersect(tmp_our, validated_genes))
# 
# length(intersect(tmp_pfc35, validated_genes))/length(tmp_pfc35)
# length(intersect(tmp_all, validated_genes))/length(tmp_all)
# length(intersect(tmp_our, validated_genes))/length(tmp_our)
# 

rm(list = c("validated_genes"))

save.image(file = paste0(save_filepath, "/step7_results.RData"))

