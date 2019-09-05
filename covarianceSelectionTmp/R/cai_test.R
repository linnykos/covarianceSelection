#' Cai Test
#'
#' @param x data matrix
#' @param y data matrix
#' @param trials number of trials
#' @param cores number of cores
#' @param prob parameter for robustness
#'
#' @return numeric (p-value)
#' @export
cai_test <- function(x, y, trials = 100, cores = 1, prob = 1){
  if(ncol(x) != ncol(y)) stop("x and y have different number of dimensions")
  if(!is.matrix(x) | !is.numeric(x)) stop("x is not a numeric matrix")
  if(!is.matrix(y) | !is.numeric(y)) stop("y is not a numeric matrix")

  doMC::registerDoMC(cores = cores)

  n1 <- nrow(x); n2 <- nrow(y)

  denom_x <- .compute_variance(x); denom_y <- .compute_variance(y)
  num_x <- .compute_sigma(x); num_y <- .compute_sigma(y)
  t_org <- .compute_covStat(num_x, num_y, denom_x, denom_y, prob = prob)

  func <- function(i){
    set.seed(i*10)
    g_x <- stats::rnorm(n1); g_y <- stats::rnorm(n2)
    boot_num_x <- .compute_bootSigma(x, g_x); boot_num_y <- .compute_bootSigma(y, g_y)
    .compute_covStat(boot_num_x, boot_num_y, denom_x, denom_y)
  }

  if(is.na(cores)){
    t_boot <- numeric(trials)
    for(i in 1:trials){ t_boot[i] <- func(i) }
  } else {
    i <- 0 #debugging purposes
    t_boot <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i)))
  }

  length(which(abs(t_boot) >= abs(t_org)))/trials
}

.compute_covStat <- function(num_x, num_y, denom_x, denom_y, squared = T){
  stopifnot(length(num_x) == length(num_y), length(denom_x) == length(denom_y))
  stopifnot(length(denom_x) == 1 | length(denom_x) == length(num_x))

  if(squared){
    res <- (num_x - num_y)^2/(denom_x + denom_y)
  } else {
    res <- abs(num_x - num_y)
  }

  max(abs(res))
}

#' Compute the bootstrap empirical covariance (numerator)
#'
#' @param mat data matrix
#' @param noise_vec vector of noise
#' @param cov_mat covariance matrix
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_bootSigma <- function(mat, noise_vec, cov_mat){
  if(length(noise_vec) != nrow(mat)) stop("length(noise_vec) not equal to nrow(mat)")

  n <- nrow(mat)
  
  mat2 <- noise_vec/n * mat
  res <- Matrix::crossprod(mat, mat2)- (sum(noise_vec)/n)*cov_mat
}

#' Compute the empirical covariance (numerator)
#'
#' @param mat data matrix
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_sigma <- function(mat){
  n <- nrow(mat)
  (n-1)/n * stats::cov(mat)
}


#' Compute the variance of the empirical covariance estimate
#'
#' @param mat data matrix
#' @param cov_mat optional argument of the empirical covariance
#' @param idx vector of the indices of the lower triangle
#'
#' @return vector
.compute_variance <- function(mat, cov_mat = NA){
  n <- nrow(mat)
  if(any(is.na(cov_mat))){
    cov_mat <- (n-1)/n*stats::cov(mat)
  } else {
    stopifnot(ncol(cov_mat) == ncol(mat))
  }
  
  mat2 <- mat^2

  crossprod(mat2)/n - cov_mat^2
}
