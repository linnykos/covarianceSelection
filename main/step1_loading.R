#format the brainspan dataset
load("../data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)

#determine brain-expressed genes
brain_expression <- covarianceSelection::brain_expression
brain_genes <- brain_expression$Gene[brain_expression$Brain_expressed != 'No']
idx <- which(colnames(genexp) %in% brain_genes)
genexp <- genexp[,idx]

#translate into synonyms
vec <- colnames(genexp)
vec <- covarianceSelection::symbol_synonyms(vec)
idx <- which(is.na(vec))
genexp <- genexp[,-idx]; vec <- vec[-idx]
colnames(genexp) <- vec

#remove duplicated genes
genexp <- genexp[,-which(duplicated(colnames(genexp)))]

#remove samples from subregions that we don't have a region for
region_subregion <- covarianceSelection::region_subregion
vec <- rownames(genexp)
subregion <- unlist(strsplit(vec,"\\."))[seq(2, length(vec)*4, 4)]
idx <- which(subregion %in% region_subregion$subregion)
genexp <- genexp[idx,]

####

#load tada dataset
tada <- covarianceSelection::tada
vec <- tada$Gene
vec <- covarianceSelection::symbol_synonyms(vec)
idx <- which(is.na(vec))
tada <- tada[-idx,]; vec <- vec[-idx]
tada$Gene <- vec

#remove duplicated tada by keeping the one with the lowest p-value
tada <- tada[-which(duplicated(tada$Gene)),]

#match the order in both datasets
idx <- which(colnames(genexp) %in% tada$Gene)
genexp <- genexp[,idx]
idx <- which(tada$Gene %in% colnames(genexp))
tada <- tada[idx,]
idx <- covarianceSelection::matching(tada$Gene, colnames(genexp))
genexp <- genexp[,idx]

dat_list <- covarianceSelection::extractor(genexp)
dat_list <- lapply(dat_list, as.matrix)

if(verbose) print(paste0("Dimension of genexp is: ", paste0(dim(genexp), collapse = ", ")))

#cleanup
rm(list = c("brain_expression", "brain_genes", "idx", "vec", "region_subregion",
            "subregion", "genexp"))

save.image(file = paste0(save_filepath, "/step1_res.RData"))
