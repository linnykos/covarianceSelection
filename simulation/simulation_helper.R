.normalize_mat <- function(mat, normalize = T){
  mat <- (mat + t(mat))/2

  eig <- eigen(mat)
  eig$values[eig$values <= 0.01] <- 0.01

  mat <- eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
  mat/max(mat)
}

.shuffle <- function(mat, percentage = 0.1, normalize = T){
  if(ncol(mat) != nrow(mat)) stop("mat is not square")
  if(percentage == 0) return(mat)

  d <- ncol(mat)
  idx <- sample(1:d, floor(percentage * d))
  mat[sort(idx), ] <- mat[idx, ]
  mat[-idx, sort(idx)] <- mat[-idx, idx]

  .normalize_mat(mat, normalize)
}

.generate_block <- function(d){
  mat <- matrix(0, d, d)
  for(j in 1:d){
    for(i in 1:j){
      mat[i,j] <- i*(d-j+1)/(d+1)
      mat[j,i] <- mat[i,j]
    }
  }

  .normalize_mat(mat)
}

.spectral_norm_mat <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  max(svd(mat1 - mat2)$d)
}
