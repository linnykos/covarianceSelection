if(verbose) print(paste0(Sys.time(), "Start of step 7: Compiling the results"))
validated_genes <- covarianceSelection::validated_genes$Gene

num_pfc35 <- length(intersect(genes_pfc35, validated_genes))  
num_nodawn <- length(intersect(genes_nodawn, validated_genes)) 
num_all <- length(intersect(genes_all, validated_genes)) 
num_our <- length(intersect(genes_our, validated_genes)) 

c(num_pfc35, num_nodawn, num_all, num_our)
c(length(genes_pfc35), length(genes_nodawn), length(genes_all), length(genes_our))
c(hmrf_pfc35$b, hmrf_all$b, hmrf_our$b)
c(hmrf_pfc35$c, hmrf_all$c, hmrf_our$c)
c(sum(as.matrix(adj_pfc35))/2, sum(as.matrix(adj_all))/2, sum(as.matrix(adj_our))/2)

cbind(names(dat_list[idx_our]), sapply(dat_list[idx_our], nrow))

# cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_pfc35, scale_vec_pfc35)[20:30,]
# cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_all, scale_vec_all)[20:30,]
# cbind(exp(seq(log(0.01), log(0.35), length.out = 30)), edges_vec_our, scale_vec_our)[20:30,]

covarianceSelection::binning(names(dat_list)[idx_our])

length(intersect(genes_pfc35, genes_our))
length(intersect(genes_pfc35, genes_nodawn))
length(intersect(genes_our, genes_nodawn))

idx_our <- which(tada$Gene %in% genes_our)
z_our <- 1 - qnorm(tada$pval.TADA[idx_our])
idx_pfc35 <- which(tada$Gene %in% genes_pfc35)
z_pfc35 <- 1 - qnorm(tada$pval.TADA[idx_pfc35])
vioplot(z_pfc35, z_our)
############

zz <- genes_our[which(!genes_our %in% genes_pfc35)]
length(zz)
length(intersect(zz, validated_genes)); zz[which(zz %in% validated_genes)]; zz[which(!zz %in% validated_genes)]
idx_our <- which(tada$Gene %in% zz)
z_our <- 1 - qnorm(tada$pval.TADA[idx_our])
quantile(z_our, probs = c(0.1,0.9))

zz <- genes_pfc35[which(!genes_pfc35 %in% genes_our)]
length(zz)
length(intersect(zz, validated_genes)); zz[which(zz %in% validated_genes)]
idx_pfc35 <- which(tada$Gene %in% zz)
z_pfc35 <- 1 - qnorm(tada$pval.TADA[idx_pfc35])
quantile(z_pfc35, probs = c(0.1,0.9))

vioplot(z_pfc35, z_our)

rm(list = c("validated_genes"))

save.image(file = paste0(save_filepath, "/step7_results", filepath_suffix, ".RData"))

