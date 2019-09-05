rm(list=ls())
library(simulation)
library(covarianceSelection)

paramMat <- cbind(10, 2, 2, 500, 10, c(0, 0.25, 0.5, 1))
colnames(paramMat) <- c("group1", "group2", "group3", "n", "d", "kappa")

# collect all the marginal densities
load("../data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)
set.seed(10)
idx <- sample(1:ncol(genexp), paramMat[1,"d"])

den_list <- lapply(idx, function(i){stats::density(genexp[,i])})

#############

generate_covariance1 <- function(vec){
  mat <- matrix(0.5, ncol = vec["d"], nrow = vec["d"])
  diag(mat) <- 1
  mat
}

generate_covariance2 <- function(vec){
  alpha <- 0.5+(0.95*vec["kappa"])*0.5
  beta <- 0.5-(0.95*vec["kappa"])*0.5
  mat <- matrix(beta, ncol = vec["d"], nrow = vec["d"])
  
  d <- vec["d"]
  d2 <- round(d/2)
  mat[1:d2, 1:d2] <- alpha
  mat[(d2+1):d, (d2+1):d] <- alpha
  diag(mat) <- 1
  
  mat
}

generate_covariance3 <- function(vec){
  generate_covariance1(vec)
}

find_cliques <- function(len, indices, threshold = 0.95){
  edges <- utils::combn(len, 2)
  
  g <- igraph::graph.empty(n = len, directed = F)
  g <- igraph::add_edges(g, edges = edges[,indices])
  
  if(igraph::ecount(g) == 0) return(NA)
  
  res <- covarianceSelection::clique_selection(g, threshold = threshold,
                                               mode = "or", verbose = F,
                                               time_limit = 300)
  
  covarianceSelection::select_clique(res, c(1:5,11,13), igraph::as_adj(g))
}

################

rule <- function(vec){
  covar1 <- generate_covariance1(vec)
  covar2 <- generate_covariance2(vec)
  covar3 <- generate_covariance3(vec)
  
  cov_list <- c(lapply(1:vec["group1"], function(x){covar1}), 
                lapply(1:vec["group2"], function(x){covar2}), 
                lapply(1:vec["group3"], function(x){covar3}))
  
  dat_list <- lapply(1:length(cov_list), function(x){
    MASS::mvrnorm(vec["n"], mu = rep(0, vec["d"]), Sigma = cov_list[[x]])
  })
  
  dat_list <- lapply(1:length(dat_list), function(x){
    covarianceSelection::nonparanormal_transformation(dat_list[[x]], den_list, 
                                                      mean_vec = rep(0, vec["n"]),
                                                      sd_vec = sqrt(diag(cov_list[[x]])))
  })
  
  for(i in (vec["group1"]+vec["group2"]+1):length(dat_list)){
    dat_list[[i]] <- (1+2*vec["kappa"])*dat_list[[i]]
  }
  
  dat_list
}

#############

y <- 10
set.seed(y)
vec <- paramMat[4,]
dat <- rule(vec)

# code modified from https://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/
scatterhist <- function(x, y, xlab="", ylab=""){
  zones <-matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  x_den <- stats::density(x)
  y_den <- stats::density(y)
  
  # manually compute what the ranges should be
  x_range <- range(x); y_range <- range(y)
  width <- max(diff(x_range), diff(y_range))
  xlim <- mean(x_range) + c(-1,1)*width/2
  ylim <- mean(y_range) + c(-1,1)*width/2
  
  par(mar=c(5,5,1,1))
  plot(x, y, pch = 16, xlim = xlim, ylim = ylim,
       col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "", ylab= "")
  par(mar=c(0,5,1,1))
  plot(NA, axes=FALSE, space=0, xlim = xlim, ylim = c(0, max(x_den$y)),
       main = "", ylab = "", xlab = "")
  polygon(c(x_den$x, x_den$x[length(x_den$x)], x_den$x[1]),
          c(x_den$y, 0, 0), col = "gray")
  par(mar=c(5,0,1,1))
  plot(NA, axes=FALSE, space=0, xlim = c(0, max(y_den$y)), ylim = ylim,
       main = "", ylab = "", xlab = "")
  polygon(c(y_den$y, 0, 0), 
          c(y_den$x, y_den$x[length(y_den$x)], y_den$x[1]),
          col = "gray")
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
  
  invisible()
}

png("../figures/scatterplot_1.png", height = 1200, width = 1200, res = 300, units = "px")
i1 <- 4; i2 <- 6
scatterhist(dat[[1]][,i1], dat[[1]][,i2], xlab = paste0("Gene ",i1), 
            ylab = paste0("Gene ", i2))
graphics.off()

png("../figures/scatterplot_2.png", height = 1200, width = 1200, res = 300, units = "px")
i1 <- 1; i2 <- 8
scatterhist(dat[[1]][,i1], dat[[1]][,i2], xlab = paste0("Gene ",i1), 
            ylab = paste0("Gene ", i2))
graphics.off()
