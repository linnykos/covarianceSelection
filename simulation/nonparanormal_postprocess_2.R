rm(list = ls())
load("../results/nonparanormal.RData")

# remove res that errored
for(i in 1:length(res)){
  rm_idx <- which(sapply(res[[i]], length) == 1)
  if(length(rm_idx) > 0) res[[i]] <- res[[i]][-rm_idx]
}

# lapply(res, function(x){ # 1 through 4
#   sapply(x, function(y){ # 1 through 25
#     sapply(y$indices_list, function(z){
#       length(z)
#     })
#   })
# })

sequence_vec <- res[[4]][[1]]$indices_list

num_mat <- sapply(1:length(sequence_vec), function(i){
  set.seed(10)

  n <- sum(paramMat[1,1:3])
  g <- igraph::graph.empty(n = n, directed = F)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::add_edges(g, edges = combn_mat[,sequence_vec[[i]]])
  
  tmp <- covarianceSelection::clique_selection(g, threshold = 0.95)
  our_num <- length(tmp[[1]])
  
  tmp <- spectral_selection(g, threshold = 0.95)
  spectral_num <- length(tmp)
  
  tmp <- covarianceSelection::tsourakakis_2013(g, threshold = 0.95)
  tsourakakis_2013_num <- length(tmp)
  
  tmp <- covarianceSelection::chen_2010(g, threshold = 0.95)
  chen_num <- length(tmp)

  c(length(sequence_vec[[i]]), our_num, spectral_num, tsourakakis_2013_num, chen_num)
})


png("../figures/figure_6.png", height = 1400, width = 2500, res = 300, units ="px")
par(mfrow = c(1,2), mar = c(5,4,4,1))

par(mfrow = c(1,2))
plot(NA, xlim = range(num_mat[1,]), ylim = range(num_mat[2:5,]), xlab = "Number of accepted hypotheses",
     ylab = "Number of selected partitions") 
#lines(num_mat[1,], num_mat[2,])
points(num_mat[1,], num_mat[2,], pch = 16, cex = 1.5)

#lines(num_mat[1,], num_mat[3,], col = rgb(205,40,54, max = 255))
points(num_mat[1,], num_mat[3,], pch = 21, col = rgb(205,40,54, max = 255), bg = "white", lwd = 2, cex = 1.5)

legend("topleft", c("Largest quasi-clique", "Spectral clustering"),
       bty="n", fill= c("black", rgb(205,40,54, max = 255)))


plot(NA, xlim = range(num_mat[1,]), ylim = range(num_mat[2:5,]), xlab = "Number of accepted hypotheses",
     ylab = "Number of selected partitions") 
points(num_mat[1,], num_mat[4,], pch = 16, col = rgb(106,164,248, max = 255), cex = 1.5)
points(num_mat[1,], num_mat[5,], pch = 21, col = rgb(149,220,144, max = 255), bg = "white", lwd = 2, cex = 1.5)

legend("topleft", c("Tsourakakis et al. (2013)", "Chen and Saad (2010)"),
       bty="n", fill= c(rgb(106,164,248, max = 255), rgb(149,220,144, max = 255)))
graphics.off()