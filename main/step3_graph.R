pthres_screening <- 0.1 # from https://arxiv.org/pdf/1506.00728.pdf, bottom of page 22
num_genes <- 6670 # from https://arxiv.org/pdf/1506.00728.pdf, bottom of page 22

#recombine dataset
selected_list <- lapply(1:length(selected_list), function(x){
  nam <- unlist(strsplit(names(selected_list)[x], "\\."))
  nam <- paste0(nam[2], ".", nam[3])
  rownames(selected_list[[x]]) <- rep(nam, nrow(selected_list[[x]]))
  selected_list[[x]]
})
dat <- do.call(rbind, selected_list)
uniq <- unique(rownames(dat))
selected_list <- lapply(uniq, function(x){
  tmp <- dat[which(rownames(dat) == x),,drop = F]
  scale(tmp, center = T, scale = F)
})
dat <- do.call(rbind, selected_list)

if(verbose) print(paste0(Sys.time(), ": Dimension of dat is: ", paste0(dim(dat), collapse = ", ")))

#screen
screened_idx <- covarianceSelection::screen(dat, tada$pval.TADA, pthres = pthres_screening,
                                       num_genes = num_genes)
dat <- dat[,unlist(screened_idx)]; tada <- tada[unlist(screened_idx),]
num_primary <- length(screened_idx$primary)
dat <- scale(dat, center = T, scale = F)

if(verbose) print(paste0(Sys.time(), ": Dimension of dat is: ", paste0(dim(dat), collapse = ", ")))
if(verbose) print(paste0("num_primary is: ", num_primary))

#store GGM results
filename_func <- covarianceSelection::filename_closure(paste0(save_filepath, "beta_vec", additional_name, "/"))
tmp <- covarianceSelection::graphicalModel_store(dat, filename_func,
                                            num_primary = num_primary, cores = cores)

#cleanup
rm(list = c("screened_idx", "tmp", "selected_list"))

save.image(file = paste0(save_filepath, "step3_res", additional_name, ".RData"))
