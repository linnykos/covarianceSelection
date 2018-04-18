# seeing if varying alpha along a single dataset is correct
# conclusion: it is

rm(list=ls())
#load("../figures/clique_simulation_RoC_2017-11-14_curve.RData")
load("../figures/clique_simulation_RoC_2017-11-23_curve.RData")
source("../figures/figure_4_helper.R")

################################

# do a preliminary check to see that the hypotheses being accepted are monotonic in alpha
bool_vec <- vector("list", 5)
for(i in 1:5){
  len <- length(curve_list[[i]])
  bool_vec[[i]] <- sapply(1:(len-1), function(x){
    all(curve_list[[i]][[1]][[x+1]] %in% curve_list[[i]][[1]][[x]])
  })
  stopifnot(all(bool_vec[[i]]))
}

###############################

combn_null <- cbind(combn(15,2), (combn(5,2)+15), (combn(5,2)+20))
idx_null <- combn_null[1,]+25*combn_null[2,]
combn_mat <- combn(25,2)
idx_all <- combn_mat[1,]+25*combn_mat[2,]
idx <- which(idx_all %in% idx_null)

hyp_tpr <- vector("list", 5)
for(i in 1:5){
  hyp_tpr[[i]] <- sapply(curve_list[[i]][[1]], function(z){
    length(which(idx %in% z))
  })
}
hyp_tpr <- lapply(hyp_tpr, function(x){x/length(idx_null)})

hyp_fpr <- vector("list", 5)
for(i in 1:5){
  hyp_fpr[[i]] <- sapply(curve_list[[i]][[1]], function(z){
    length(which(!z %in% idx))
  })
}

hyp_fpr <- lapply(hyp_fpr, function(x){x/(length(idx_all) - length(idx_null))})

##################

hyp_tpr <- lapply(hyp_tpr, function(x){
  c(x, 0)
})
hyp_fpr <- lapply(hyp_fpr, function(x){
  c(x, 0)
})

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:5){
  lines(hyp_fpr[[i]], hyp_tpr[[i]], lwd = 2, col = i)
}

lines(c(0,1), c(0,1), lty = 2, lwd = 2, col = "red")

###################################
#try again for the cliques

indices <- vector("list", 5)
edges <- combn(25,2)

for(i in 1:5){
  set.seed(10)
  indices[[i]] <- lapply(curve_list[[i]][[1]], function(x){
    largest_clique(25, edges[,x], c(1:3,16,21), 0.95)
  })
}

idx_tpr <- vector("list", 5)
for(i in 1:5){
  idx_tpr[[i]] <- sapply(indices[[i]], function(x){
    length(which(1:15 %in% x))/15
  })
}

idx_fpr <- vector("list", 5)
for(i in 1:5){
  idx_fpr[[i]] <- sapply(indices[[i]], function(x){
    length(which(!x %in% 1:15))/(25-15)
  })
}

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
for(i in 1:5){
  lines(idx_fpr[[i]], idx_tpr[[i]], lwd = 2, col = i)
}

lines(c(0,1), c(0,1), lty = 2, lwd = 2, col = "red")


#######################################

## zoom in some particular settings
edges1 <- combn(25, 2)[,curve_list[[4]][[1]][[16]]]
edges2 <- combn(25, 2)[,curve_list[[4]][[1]][[17]]]


#compute the adjacency matrix
adj1 <- matrix(0, 25, 25)
diag(adj1) <- 1
for(i in 1:ncol(edges1)){
  adj1[edges1[1,i], edges1[2,i]] <- 1
  adj1[edges1[2,i], edges1[1,i]] <- 1
}

adj2 <- matrix(0, 25, 25)
diag(adj2) <- 1
for(i in 1:ncol(edges2)){
  adj2[edges2[1,i], edges2[2,i]] <- 1
  adj2[edges2[2,i], edges2[1,i]] <- 1
}

.clockwise90 = function(a) { t(a[nrow(a):1,]) }

#colors
pale <- rgb(247, 234, 200, max = 255)
red <- rgb(219, 51, 64, max = 255)

par(mfrow = c(1,2), mar = c(3,3,2,1))
image(.clockwise90(adj1), asp = T, xlab = "Index locations", ylab = "Index locations",
      main = "Adjacency matrix",
      xaxt = "n", bty = "n", yaxt = "n", col = c(pale, red), mgp = c(1,0,0))
image(.clockwise90(adj2), asp = T, xlab = "Index locations", ylab = "Index locations",
      main = "Adjacency matrix",
      xaxt = "n", bty = "n", yaxt = "n", col = c(pale, red), mgp = c(1,0,0))
