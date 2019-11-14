fdr_cutoff_alt <- 0.1

if(verbose) print(paste0(Sys.time(), "Start of step 2: No DAWN analysis"))

genes_nodawn <- sort(as.character(tada[which(tada$qvalue <= fdr_cutoff_alt),"Gene"]))

save.image(file = paste0(save_filepath, "/step2_nodawn_analysis", filepath_suffix, ".RData"))