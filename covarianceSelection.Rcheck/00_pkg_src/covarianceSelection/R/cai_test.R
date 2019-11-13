#' Cai Test
#'
#' @param x data matrix
#' @param y data matrix
#' @param trials number of trials
#' @param prob numeric between 0 and 1
#'
#' @return numeric (p-value)
#' @export
cai_test <- function(x, y, trials = 100, prob = 1){
  if(ncol(x) != ncol(y)) stop("x and y have different number of dimensions")
  if(!is.matrix(x) | !is.numeric(x)) stop("x is not a numeric matrix")
  if(!is.matrix(y) | !is.numeric(y)) stop("y is not a numeric matrix")

  diag_idx <- which(lower.tri(diag(ncol(x)), diag = T))

  n1 <- nrow(x); n2 <- nrow(y)

  num_x <- .compute_sigma(x, diag_idx); num_y <- .compute_sigma(y, diag_idx)
  denom_x <- .compute_variance(x, num_x, diag_idx); denom_y <- .compute_variance(y, num_y, diag_idx)
  t_org <- .compute_covStat(num_x, num_y, denom_x, denom_y, prob = prob)

  func <- function(i){
    set.seed(i*10)
    g_x <- stats::rnorm(n1); g_y <- stats::rnorm(n2)
    boot_num_x <- .compute_bootSigma(x, g_x, num_x, diag_idx)
    boot_num_y <- .compute_bootSigma(y, g_y, num_y, diag_idx)
    .compute_covStat(boot_num_x, boot_num_y, denom_x, denom_y, prob = 1)
  }

  i <- 0 #debugging purposes
  t_boot <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i)))

  length(which(abs(t_boot) >= abs(t_org)))/trials
}

.compute_covStat <- function(num_x, num_y, denom_x, denom_y, squared = T, prob = 1){
  if(squared){
    res <- (num_x - num_y)^2/(denom_x + denom_y)
  } else {
    res <- abs(num_x - num_y)
  }

  if(prob == 1) max(abs(res)) else {
    stats::quantile(abs(res), prob = prob)
  }
}

#' Compute the bootstrap empirical covariance (numerator)
#'
#' @param mat data matrix
#' @param noise_vec vector of noise
#' @param cov_vec vectorized covariance matrix of length \code{idx}
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_bootSigma <- function(mat, noise_vec, cov_vec, idx){
  n <- nrow(mat)
  
  mat2 <- noise_vec/n * mat
  Matrix::crossprod(mat, mat2)[idx] - (sum(noise_vec)/n)*cov_vec
}

#' Compute the empirical covariance (numerator)
#'
#' @param mat data matrix
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_sigma <- function(mat, idx){
  n <- nrow(mat)
  (n-1)/n * stats::cov(mat)[idx]
}


#' Compute the variance of the empirical covariance estimate
#'
#' @param mat data matrix
#' @param cov_vec vectorized covariance matrix of length \code{idx}
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_variance <- function(mat, cov_vec, idx){
  n <- nrow(mat)
  mat2 <- mat^2

  crossprod(mat2)[idx]/n - cov_vec^2
}
