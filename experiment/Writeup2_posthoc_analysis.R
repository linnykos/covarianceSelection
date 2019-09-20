rm(list=ls())
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection.RData")
prob = 0.9999

dat_list <- lapply(dat_list, scale, center = T, scale = F)
diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
len <- length(dat_list)
combn_mat <- utils::combn(len, 2)

num_list <- lapply(dat_list, function(x){covarianceSelection:::.compute_sigma(x, diag_idx)})
denom_list <- covarianceSelection:::.compute_all_denom(dat_list, num_list, diag_idx)

ncores <- 20
doMC::registerDoMC(cores = ncores)

func <- function(x){
  print(x)
  covarianceSelection:::.compute_covStat(num_list[[combn_mat[1,x]]], num_list[[combn_mat[2,x]]],
                                         denom_list[[combn_mat[1,x]]], denom_list[[combn_mat[2,x]]],
                                         squared = T, prob = prob)
}

t_vec <- foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), func(i))

save(t_vec, file = "/raid6/Kevin/covarianceSelection/results/step3_subjectselection_updated.RData")

################################

rm(list=ls())
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection.RData")
load("/raid6/Kevin/covarianceSelection/results/step3_subjectselection_updated.RData")
t_vec <- unlist(t_vec)

stepdown_res <- lapply(seq(0, 1, length.out = 21), function(alpha){
  print(alpha)
  return_pvalue = F
  verbose = F
  stopifnot(class(stepdown_obj) == "stepdown")
  len <- length(t_vec)
  idx_all <- rep(TRUE, len)
  
  round <- 1
  while(TRUE){
    if(sum(idx_all) == 0) break()
    
    boot_mat <- stepdown_obj$boot[[round]][, which(idx_all), drop = F]
    t_boot <- apply(boot_mat, 1, function(x){abs(max(x))})
    
    cutoff <- stats::quantile(abs(t_boot), 1-alpha)
    idx <- intersect(which(abs(t_vec) >= cutoff), which(idx_all))
    
    if(length(idx) == 0) break()
    
    idx_all[idx] <- FALSE
    if(verbose) print(paste0("In stepdown, finished round ", round, " with ",
                             sum(idx_all), " null hypothesis remaining"))
    round <- round+1
    
    if(round > length(stepdown_obj$boot)) {
      warning("Function ran out of bootstrap rounds to work with.")
      break()
    }
  }
  
  which(idx_all)
})

##########################

idx <- stepdown_res[[3]]
# let's construct a graph
n <- length(dat_list)
g <- igraph::graph.empty(n = n, directed = F)
combn_mat <- utils::combn(length(dat_list), 2)
g <- igraph::add_edges(g, edges = combn_mat[,idx])
save.image("/raid6/Kevin/covarianceSelection/results/step3_subjectselection_updated2.RData")

