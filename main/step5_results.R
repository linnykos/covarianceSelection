iossifov <- covarianceSelection::iossifov
iossifov <- covarianceSelection::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

fisher <- covarianceSelection::enrichment_test(autism_genes, iossifov, report$Gene)
fisher$percentage <- fisher$contigency[1,1]/(fisher$contigency[1,1]+fisher$contigency[1,2])

scale_free <- covarianceSelection::compute_scale_free(adj_gene)

g_gene <- covarianceSelection::adj_to_graph(adj_gene)

higher_crit <- fdrtool::hc.thresh(diagnostic_vec, alpha0 = 0.5)

summary_results <- list(percentage = fisher$percentage,
                        num_aut = length(autism_genes), bin = bin,
                        num_partition = length(selected_names),
                        num_samples = nrow(dat), num_genes = ncol(dat),
                        scale_free = scale_free, hmrf_c = hmrf_res$c,
                        higher_crit = higher_crit,
                        diagnostic = diagnostic_vec)

rm(list = c("iossifov", "bin", "cores", "scale_free", "dat", "numerator",
            "denominator", "clustering_dist", "tada"))

save.image(file = paste0(save_filepath, "step5_res", additional_name, ".RData"))
