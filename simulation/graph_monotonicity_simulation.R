rm(list=ls())
library(simulation)
library(covarianceSelection)

paramMat <- as.matrix(expand.grid(50, 100, c(1:4)))
colnames(paramMat) <- c("n", "spacing", "method")

################

rule <- function(vec){
  invisible()
}

criterion <- function(dat, vec, y){
  combn_mat <- utils::combn(vec["n"], 2)
  g <- igraph::graph.empty(n = vec["n"], directed = F)
  g <- igraph::add_edges(g, edges = combn_mat)
  num_edges <- ncol(combn_mat)
  remove_num <- floor(num_edges / (vec["spacing"]+1))
  
  i <- 1
  res <- rep(NA, vec["spacing"])
  while(i <= vec["spacing"]){
    print(i)
    if(vec["method"] == 1){
      res[i] <- length(covarianceSelection::tsourakakis_2013(g))
    } else if(vec["method"] == 2){
      res[i] <- length(covarianceSelection::chen_2010(g))
    } else if(vec["method"] == 3){
      res[i] <- length(covarianceSelection::anderson_2009(g))
    } else {
      res[i] <- length(covarianceSelection::tsourakakis_2014_approximate(g))
    }
    
    # remove indices for the next iteration
    idx <- sample(1:igraph::ecount(g), remove_num)
    g <- igraph::delete_edges(g, idx)
    
    i <- i+1
  }
}