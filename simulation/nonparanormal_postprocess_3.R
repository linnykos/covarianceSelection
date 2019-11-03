rm(list = ls())
load("../results/nonparanormal.RData")

# remove res that errored
for(i in 1:length(res)){
  rm_idx <- which(sapply(res[[i]], length) == 1)
  if(length(rm_idx) > 0) res[[i]] <- res[[i]][-rm_idx]
}

set.seed(10)
ncores <- 21
doMC::registerDoMC(cores = ncores)

func <- function(z){
  set.seed(10)
  
  n <- sum(paramMat[1,1:3])
  g <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::add_edges(g, edges = combn_mat[,z])
  tmp <- covarianceSelection::clique_selection(g, threshold = 0.95, verbose = F, time_limit = 60)
  
  if(length(res) > 1){
    len <- sapply(tmp, function(zz){length(intersect(zz, 1:paramMat[1,1]))})
    tmp <- tmp[[which.max(len)]]
  } else {
    tmp <- tmp[[1]]
  }
  tmp
}


partition_mat <- lapply(res, function(x){
  cat("\n")
  lapply(x, function(y){
    cat("*")
    foreach::"%dopar%"(foreach::foreach(i = 1:length(y$indices_list)), func(y$indices_list[[i]]))
  })
})

save.image("../results/nonparanormal_3.RData")
