rm(list=ls())
library(longitudinalGM)
alpha_tuning_seq <- seq(0, 0.35, by = 0.025)

iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

sequence_res <- vector("list", length(alpha_tuning_seq))

load("../results/step5_res_pfc35.RData")
gene_intersect <- intersect(iossifov, autism_genes)
pfc_res <- list(fisher = fisher$contigency, bin = summary_results$bin,
                num_samples = summary_results$num_samples,
                gene_intersect = gene_intersect)

for(idx in 1:length(alpha_tuning_seq)){
  load(paste0("../results/step5_res_alpha", alpha_tuning_seq[idx], ".RData"))

  gene_intersect <- intersect(iossifov, autism_genes)

  sequence_res[[idx]] <- list(fisher = fisher$contigency, bin = summary_results$bin,
                              num_samples = summary_results$num_samples,
                              autism_genes = autism_genes,
                              gene_intersect = gene_intersect)
}
names(sequence_res) <- paste0("alpha_", alpha_tuning_seq)

############
n <- sum(pfc_res$fisher[1,])
p <- pfc_res$fisher[1,1]/n
std <- sqrt(n*p*(1-p))

#############

# find the union
start <- 3
end <- 15
zz <- sequence_res[[start]]$gene_intersect
for(i in start:end){
  zz <- c(zz, sequence_res[[i]]$gene_intersect)
}
zz <- sort(unique(zz))

count_vec <- rep(0, length(zz))
names(count_vec) <- zz

for(i in start:end){
  idx <- which(names(count_vec) %in% sequence_res[[i]]$gene_intersect)
  count_vec[idx] <- count_vec[idx] + 1
}

length(which(count_vec >= 8))
sort(names(which(count_vec >= 8)))

###########

# find the union
start <- 3
end <- 15
zz <- as.character(sequence_res[[start]]$autism_genes)
for(i in start:end){
  zz <- c(zz, as.character(sequence_res[[i]]$autism_genes))
}
zz <- sort(unique(zz))

count_vec <- rep(0, length(zz))
names(count_vec) <- zz

for(i in start:end){
  idx <- which(names(count_vec) %in% as.character(sequence_res[[i]]$autism_genes))
  count_vec[idx] <- count_vec[idx] + 1
}

length(which(count_vec >= 8))
sort(names(which(count_vec >= 8)))

