trials <- 200
num_permutations <- 500
quant <- 1
num_genes_partition <- 200
threshold <- 0.95

idx <- which(sapply(dat_list, function(x){ifelse(nrow(x) >= 5, T, F)}))
dat_list <- dat_list[idx]
tmp_list <- lapply(dat_list, function(x){x[,1:min(ncol(x), num_genes_partition),drop=F]})

#select pfc-35
len <- length(dat_list)

#in-group tests
if(additional_name == "" | additional_name == "_nodenom"){
  if(verbose) print(paste0(Sys.time(), ": Starting to use stepdown"))

  selected_p <- covarianceSelection::stepdown(tmp_list, trials = trials,
                                         denominator = ifelse(additional_name == "", T, F),
                                         cores = cores,
                                         alpha = alpha, verbose = verbose)

  g <- igraph::graph.empty(n = len, directed = F)
  g <- igraph::add_edges(g, edges = combn(len, 2)[,selected_p])
  initial_idx <- grep("PFC\\.[3-5]", names(dat_list))

  clique_list <- covarianceSelection::clique_selection(g, threshold = threshold,
                                                  target_idx = initial_idx, verbose = F)

  selected_idx <- covarianceSelection::select_clique(clique_list, initial_idx,
                                                as.matrix(igraph::as_adj(g)))

} else if (additional_name == "_alternative") {
  selected_idx <- 1:length(dat_list)
} else if (additional_name == "_pfc35"){
  selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
} else {
  stop("variable 'additional_name' was not set.")
}

if(verbose) print(paste0(Sys.time(), ": Number of selected partitions: ", length(selected_idx)))

selected_list <- dat_list[selected_idx]
num_samples_partition <- sapply(dat_list, nrow)
tmp_list <- tmp_list[selected_idx]
selected_names <- names(selected_list)

bin <- covarianceSelection::binning(selected_names)

#goodness of fit diagnostic
set.seed(10)
diagnostic_vec <- covarianceSelection::goodness_of_fit(tmp_list, num_permutations, trials = trials,
                                                  cores = cores, verbose = verbose)

if(verbose) print(paste0(Sys.time(), ": Overview of diagnostic: ",
                         paste0(round(stats::quantile(diagnostic_vec),2), collapse = ", ")))

#cleanup
rm(list = c("genexp", "dat_list", "idx", "tmp_list"))

save.image(file = paste0(save_filepath, "step2_res", additional_name, ".RData"))

