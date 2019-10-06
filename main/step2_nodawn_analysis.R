if(verbose) print(paste0(Sys.time(), "Start of step 2: No DAWN analysis"))

cutoff <- quantile(tada$qvalue, prob = 200/length(tada$qvalue), na.rm = T)
genes_nodawn <- sort(as.character(tada[which(tada$qvalue <= cutoff),"Gene"]))

rm(list = c("cutoff"))

save.image(file = paste0(save_filepath, "/step2_nodawn_analysis.RData"))