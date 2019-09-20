rm(list=ls())
tada <- read.csv("../../raw_data/TADA_Results_2231trios_1333trans_1601cases_5397controls_March26_pvalues.csv") # 18735 genes
tada <- tada[,which(colnames(tada) %in% c("Gene", "dn.LoF", "qvalue", "pval.TADA"))]
vec <- covarianceSelection::symbol_synonyms( tada$Gene, verbose = T)
tada2 <- tada

unknown_genes_idx <- which(sapply(vec, length) == 0)
vec2 <- vec[-unknown_genes_idx]; vec2 <- unlist(vec2)

# about to do something janky
load("../results/step2_pfc35_analysis.RData")

vec_mapping <- rep(NA, length(vec))
for(i in 1:length(vec_mapping)){
  if(i %% floor(length(vec_mapping)/10) == 0) cat('*')
  if(length(vec[[i]]) > 0){
    res <- which(tada$Gene == vec[[i]])
    if(length(res) > 0){
      vec_mapping[i] <- res
    }
  }
}

q_vec <- rep(NA, max(vec_mapping, na.rm = T))
for(i in 1:length(vec_mapping)){
  if(!is.na(vec_mapping[i])){
    q_vec[vec_mapping[i]] <- tada2[i,"qvalue"]
  }
} #something weird going on... but let's go with it

cutoff <- quantile(q_vec, prob = 200/length(q_vec), na.rm = T)
genes_nodawn <- tada$Gene[which(q_vec <= cutoff)]
