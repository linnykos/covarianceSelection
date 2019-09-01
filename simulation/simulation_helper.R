.clean_mat <- function(mat, normalize = T){
  stopifnot(nrow(mat) == ncol(mat))
  # symmetrize
  mat <- (mat + t(mat))/2

  # ensure PSD
  if(normalize){
    eig <- eigen(mat)
    eig$values[eig$values <= 0.01] <- 0.01
    mat <- eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
  }
  
  mat
}

.generate_block <- function(d, alpha = 0.75, beta = .25, spillover_percentage = 0,
                            normalize = T){
  mat <- matrix(beta, d, d)
  cutoff <- round(d/2)
  
  cluster_lab <- rep(2, d); cluster_lab[1:cutoff] <- 1
  # account for spillover
  num_spillover <- round(d*spillover_percentage)
  if(num_spillover > 0){
    for(i in 1:2){
      idx <- which(cluster_lab == i)[1:num_spillover]
      cluster_lab[idx] <- 3
    }
  }
  
  for(i in 1:3){
    idx <- which(cluster_lab == i)
    if(length(idx) > 0){
      mat[idx, idx] <- alpha
    }
  }
  
  diag(mat) <- 1
  
  .clean_mat(mat, normalize = normalize)
}