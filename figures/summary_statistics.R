rm(list=ls())

statements <- vector("list")

load("../results/step1_res.RData")
statements[[1]] <- paste0("There are ", ncol(dat_list[[1]]), " genes and ",
                          sum(sapply(dat_list, nrow)), " microarray samples across ",
                          length(dat_list), " partitions, ",
                          length(which(sapply(dat_list, function(x){ifelse(nrow(x) >= 5, T, F)}))),
                          " of which have more than 5 samples.")

load("../results/step5_res.RData")
iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

######

statements[[2]] <- paste0("The parameters used in this analysis are: ",
                          "alpha = ", alpha, ", number of bootstrap trails = ",
                          trials, ", number of divisions tried in diagnostic = ",
                          num_permutations, ", number of genes in partition selection = ",
                          num_genes_partition, ", number of genes in downstream analysis = ",
                          num_genes, ", number of risk genes selected = ",
                          num_target, ".")
statements[[3]] <- paste0("There are ", length(iossifov), " genes in the independent study.")
statements[[4]] <- paste0("A total of ", sum(summary_results$bin), " of partitions, collectively ",
                          "containing ", summary_results$num_samples, " samples, were selected:")
statements[[5]] <- summary_results$bin

#########

overlap <- sort(iossifov[which(iossifov %in% as.character(autism_genes))])

statements[[6]] <- paste0("A total of ", summary_results$num_aut, " of risk genes were found ",
                          "using our current method, ",
             length(overlap), " of which are also in the independent study (",
             round(length(overlap)/summary_results$num_aut*100, 2), "%):")
statements[[7]] <- overlap

load("../results/step5_res_pfc35.RData")

overlap_pfc <- sort(iossifov[which(iossifov %in% as.character(autism_genes))])
statements[[8]] <- paste0("The previous method found ", summary_results$num_aut, " risk genes, ",
                          length(overlap_pfc), " of which were also the independent study (",
                          round(length(overlap_pfc)/summary_results$num_aut*100, 2), "%).")

#######

autism_genes_main <- autism_genes
overlap_main <- overlap

alpha_tuning_seq <- seq(0.05, 0.35, by = 0.025)
autism_list <- vector("list", length(alpha_tuning_seq))
overlap_list <- vector("list", length(alpha_tuning_seq))

for(idx in 1:length(alpha_tuning_seq)){
  load(paste0("../results/step5_res_alpha", alpha_tuning_seq[idx], ".RData"))

  overlap_list[[idx]] <- sort(iossifov[which(iossifov %in% as.character(autism_genes))])
  autism_list[[idx]] <- autism_genes
}

overlap_table <- table(unlist(overlap_list))

at_least_cutoff <- 7

statements[[9]] <- paste0("Of the ", length(overlap_main), " overlap risk genes in our analysis, ",
                          length(which(overlap_table[names(overlap_table) %in% overlap_main] >= 6)),
                          " of them appear at least ", at_least_cutoff, " of the ", length(alpha_tuning_seq), " analyses with ",
                          " alpha varying from ", min(alpha_tuning_seq), " to ", max(alpha_tuning_seq), ".")

autism_table <- table(unlist(autism_list))

statements[[10]] <- paste0("Similarly, of the ", length(autism_genes_main), " risk genes ",
                           " in our analysis, ",
                           length(which(autism_table[names(autism_table) %in% autism_genes_main] >= 6)),
                           " of them appear at least ", at_least_cutoff, " of the ", length(alpha_tuning_seq), " analyses.")

print(statements)


