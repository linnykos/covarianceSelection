if(verbose) print("Start of step 0: Loading")

#format the brainspan dataset
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp) # 1340 x 16947

#determine brain-expressed genes
brain_expression <- covarianceSelection::brain_expression
brain_genes <- brain_expression$Gene[brain_expression$Brain_expressed != 'No']
idx <- which(colnames(genexp) %in% brain_genes)
genexp <- genexp[,idx] # 1340 x 14370

#translate into synonyms
vec <- covarianceSelection::symbol_synonyms(colnames(genexp), verbose = T)
unknown_genes_idx <- which(is.na(vec))
vec <- vec[-unknown_genes_idx]
genexp <- genexp[-unknown_genes_idx] # 1340 x 14297
colnames(genexp) <- vec

#average non-unique genes
genexp <- covarianceSelection:::average_same_columns(genexp) # 1340 x 14246

#remove samples from subregions that we don't have a region for
region_subregion <- covarianceSelection::region_subregion
vec <- rownames(genexp)
subregion <- unlist(strsplit(vec,"\\."))[seq(2, length(vec)*4, 4)]
idx <- which(subregion %in% region_subregion$subregion)
genexp <- genexp[idx,] # 1294 x 14249

####

#load tada dataset
tada <- covarianceSelection::tada

#match the order in both datasets
idx <- which(colnames(genexp) %in% tada$Gene)
genexp <- genexp[,idx] # 1294 x 13964
idx <- which(tada$Gene %in% colnames(genexp))
tada <- tada[idx,] # 13964 genes
idx <- covarianceSelection::matching(tada$Gene, colnames(genexp))
genexp <- genexp[,idx]

dat_list <- covarianceSelection::extractor(genexp) # 212 partitions
dat_list <- lapply(dat_list, as.matrix, drop = F)

# remove partitions too small
idx <- which(sapply(dat_list, function(x){ifelse(nrow(x) >= 5, T, F)}))
dat_list <- dat_list[idx] # 125 partitions
dat_list <- lapply(dat_list, scale, scale = F)

if(verbose) print(paste0("Dimension of genexp is: ", paste0(dim(genexp), collapse = ", ")))

#cleanup
rm(list = c("brain_expression", "brain_genes", "idx", "vec", "region_subregion",
            "subregion", "genexp", "unknown_genes_idx"))

save.image(file = paste0(save_filepath, "/step0_loading", filepath_suffix, ".RData"))
