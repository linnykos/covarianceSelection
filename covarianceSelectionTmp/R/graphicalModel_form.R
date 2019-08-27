#' Graphical model estimate
#'
#' After a call of \code{graphicalModel_store}, \code{graphicalmodel_form}
#' will take the saved glmnet objects (dictated by \code{filename_func})
#' to form the precision matrix
#'
#' @param dat the nxd matrix
#' @param filename_func a function to create the file name. it takes in
#' \code{i}, the dimension index, as a parameter
#' @param lambda lambda value
#' @param num_primary integer from 1 to \code{ncol(dat)} that denotes
#' the number of genes to consider
#' @param tol numeric
#' @param cores number of cores
#' @param verbose boolean for verbose
#'
#' @return a precision matrix estimate, dxd
#' @export
graphicalModel_form <- function(dat, filename_func, lambda = NA, num_primary = ncol(dat),
                                tol = 1e-6, cores = 1, verbose = F){
  stopifnot(num_primary >= 1, num_primary <= ncol(dat), num_primary %% 1 == 0)
  d <- ncol(dat)

  if(verbose) print(paste0(Sys.time(), ": Extracting coefficients"))
  coef_list <- .extract_beta(filename_func, lambda, d = d, num_primary = num_primary,
                             cores = cores)
  coef_mat <- do.call(cbind, coef_list)

  if(num_primary < d) coef_mat <- .readjust_coef(coef_mat, d)
  if(verbose) print(paste0(Sys.time(), ": Computing sigma"))
  sigma_vec <- .compute_sigma_vec(dat, coef_mat, cores = cores)


  if(verbose) print(paste0(Sys.time(), ": Forming precision matrix"))
  i <- 1 #debugging reasons
  doMC::registerDoMC(cores = cores)
  func <- function(x){
    vec <- numeric(d)
    vec <- -coef_mat[,x]/sigma_vec[x]
    vec[x] <- 1/sigma_vec[x]

    vec
  }
  prec_mat <- foreach::"%dopar%"(foreach::foreach(i = 1:d), func(i))
  prec_mat <- do.call(cbind, prec_mat)

  if(verbose) print(paste0(Sys.time(), ": Symmetrizing matrix"))
  prec_mat <- .symmetrize(prec_mat)
  prec_mat[abs(prec_mat) < tol] <- 0

  prec_mat
}

.readjust_coef <- function(mat, d){
  num_primary <- ncol(mat)
  mat <- cbind(mat, matrix(0, nrow = nrow(mat), ncol = d - num_primary))
  mat[1:num_primary, (num_primary:d)] <- t(mat[(num_primary:d), 1:num_primary])
  mat
}

#' Loads the \code{glmnet} or \code{cv.glmnet} object and extracts the relevant
#' coefficients.
#'
#' The object loaded is assumed to be called \code{res}. If \code{res} is
#' not of class \code{cv.glmnet}, \code{lambda} is used. Otherwise,
#' \code{lambda} is set to \code{res$lambda.1se}.
#'
#' @param filename_func a function to create the file name. it takes in
#' \code{i}, the dimension index, as a parameter
#' @param lambda lambda value
#' @param d number of dimensions
#' @param num_primary integer from 1 to \code{ncol(dat)} that denotes
#' the number of genes to consider
#' @param cores number of cores
#'
#' @return list of vectors
.extract_beta <- function(filename_func, lambda = NA, d, num_primary = d, cores = 1){
  doMC::registerDoMC(cores = cores)
  filename_vec <- sapply(1:d, filename_func)

  res <- NA; i <- 1 #debugging reasons
  func <- function(x){
    load(filename_vec[x])
    vec <- rep(0, d)

    if(inherits(res,"cv.glmnet")) lambda <- res$lambda.1se
    stopifnot(!is.na(lambda))

    vec[-x] <- as.numeric(stats::coef(res, s = lambda))[-1]

    vec
  }

  foreach::"%dopar%"(foreach::foreach(i = 1:num_primary), func(i))
}
