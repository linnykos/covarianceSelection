# in this test, we see if aggregating the datasets by bins helps (effectively removing one entire layer of complexity)
rm(list=ls())

ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
#install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)

verbose <- T
save_filepath <- "/raid6/Kevin/covarianceSelection/results"

load("../results/step1_screening.RData") # we want the "all_genes" variable

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
unknown_genes_idx <- which(sapply(vec, length) == 0)
vec <- vec[-unknown_genes_idx]; vec <- unlist(vec)
genexp <- genexp[-unknown_genes_idx] # 1340 x 14297
colnames(genexp) <- vec

#average non-unique genes
genexp <- covarianceSelection:::average_same_columns(genexp) # 1340 x 14249 #EDIT THIS

#remove samples from subregions that we don't have a region for
region_subregion <- covarianceSelection::region_subregion
vec <- rownames(genexp)
subregion <- unlist(strsplit(vec,"\\."))[seq(2, length(vec)*4, 4)]
idx <- which(subregion %in% region_subregion$subregion)
genexp <- genexp[idx,] # 1294 x 14249

####

#load tada dataset
tada <- read.csv("../../raw_data/TADA_Results_2231trios_1333trans_1601cases_5397controls_March26_pvalues.csv") # 18735 genes
tada <- tada[,c(1,3,17)]
vec <- covarianceSelection::symbol_synonyms( tada$Gene, verbose = T)
unknown_genes_idx <- which(sapply(vec, length) == 0)
tada <- tada[-unknown_genes_idx,] # 18700 genes
vec <- vec[-unknown_genes_idx]; vec <- unlist(vec)
tada$Gene <- vec

#remove duplicated tada by keeping the one with the lowest p-value
tada <- tada[-which(duplicated(tada$Gene)),] #18498 genes

#match the order in both datasets
idx <- which(colnames(genexp) %in% tada$Gene)
genexp <- genexp[,idx] # 1294 x 13964
idx <- which(tada$Gene %in% colnames(genexp))
tada <- tada[idx,] # 13964 genes
idx <- covarianceSelection::matching(tada$Gene, colnames(genexp))
genexp <- genexp[,idx]

dat_list <- covarianceSelection::extractor(genexp) # 212 partitions
dat_list <- lapply(dat_list, as.matrix, drop = F)

# REFORMULATE THE DATA HERE
name_mat <- sapply(names(dat_list), covarianceSelection:::.split_name)
time <- covarianceSelection:::.time_to_splits(as.numeric(name_mat[3,]), splits = c(2,5,8,15))
region <- factor(as.character(name_mat[2,]), levels = c("PFC", "VIIAS", "SHA", "MDCBC"))

dat_aggregate <- vector("list", length(levels(time))*length(levels(region)))
for(i in 1:length(levels(time))){
  for(j in 1:length(levels(region))){
    idx <- intersect(which(time == levels(time)[i]), which(region == levels(region)[j]))
    dat_aggregate[[(i-1)*length(levels(region))+j]] <- do.call(rbind, dat_list[idx])
    names(dat_aggregate)[(i-1)*length(levels(region))+j] <- paste0( levels(region)[j], "-", levels(time)[i])
  }
}

idx <- which(sapply(dat_aggregate, function(x){ifelse(is.matrix(x), nrow(x), 0)}) >= 15)
dat_aggregate <- dat_aggregate[idx]

#cleanup
rm(list = c("brain_expression", "brain_genes", "idx", "vec", "region_subregion",
            "subregion", "genexp", "unknown_genes_idx", "dat_list"))

#######################

for(i in 1:length(dat_aggregate)){
  dat_aggregate[[i]] <- scale(dat_aggregate[[i]][,all_genes], center = T, scale = T)
}

trials <- 100
ncores <- 15

save(trials, file = paste0(save_filepath, "/test.RData"))
stepdown_obj <- covarianceSelection::stepdown_path(dat_aggregate, trials = trials, cores = ncores, verbose = verbose,
                                                   iterations = 7, file = paste0(save_filepath, "/step3_subjectselection_tmp2_2.RData"),
                                                   prob = 0.9)
save.image(file = paste0(save_filepath, "/step3_subjectselection_tmp_2.RData"))
stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  covarianceSelection::stepdown_choose(stepdown_obj, alpha = alpha, return_pvalue = T)
})

save.image(file = paste0(save_filepath, "/step3_subjectselection_experiment_2.RData"))

