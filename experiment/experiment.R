rm(list=ls())
load("../results/step4_ouranalysis_tmp.RData")

res_list <- vector("list", 4)
res_list[[1]] <- tsourakakis_2013(g)
res_list[[2]]  <- chen_2010(g)
res_list[[3]]  <- anderson_2009(g)
res_list[[4]]  <- tsourakakis_2014_approximate(g)

n <- igraph::vcount(g)
igraph::V(g)$name <- 1:n
for(i in 1:4){
  print(paste0("Size : ", length(res_list[[i]]), " / Density : ", 
               round(.chen_check_density(g, as.character(res_list[[i]])), 5)))
}
