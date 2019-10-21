rm(list=ls())
library(simulation)
library(covarianceSelection)
source("../simulation/graph_simulation_helper.R")

paramMat <- as.matrix(expand.grid(50, 100, c(1:4), c(1,2)))
colnames(paramMat) <- c("n", "spacing", "quasi_method", "remove_method")

################

rule <- function(vec){
  NA
}

criterion <- function(dat, vec, y){
  combn_mat <- utils::combn(vec["n"], 2)
  g <- igraph::graph.empty(n = vec["n"], directed = F)
  g <- igraph::add_edges(g, edges = combn_mat)
  num_edges <- ncol(combn_mat)
  remove_num <- floor(num_edges / (vec["spacing"]+1))
  
  i <- 1
  res <- rep(NA, vec["spacing"])
  set.seed(10*y)
  while(i <= vec["spacing"]){
    print(i)
    if(vec["quasi_method"] == 1){
      res[i] <- length(covarianceSelection::tsourakakis_2013(g))
    } else if(vec["quasi_method"] == 2){
      res[i] <- length(covarianceSelection::chen_2010(g))
    } else if(vec["quasi_method"] == 3){
      res[i] <- length(covarianceSelection::anderson_2009(g))
    } else {
      res[i] <- length(covarianceSelection::tsourakakis_2014_approximate(g))
    }
    
    # remove indices for the next iteration
    if(vec["remove_method"] == 1){
      g <- .remove_edges_uniform(g, remove_num)
    } else {
      g <- .remove_edges_proportional(g, remove_num)
    }
    
    i <- i+1
  }
}
