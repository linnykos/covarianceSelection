#' Graphical model estimate
#'
#' Estimated using neighbhorhood selection, cross validation to select lambda
#'
#' @param dat the nxd matrix
#' @param verbose boolean
#'
#' @return a precision matrix estimate, dxd
#' @export
graphicalModel <- function(dat, lambda = "lambda.1se", verbose = F, tol = 1e-6){
  n <- nrow(dat); d <- ncol(dat)

  if(verbose) print("Starting to estimate coefficients")
  coef_list <- .compute_reg_coefficients_cv(dat, lambda = lambda, verbose = verbose)
  coef_mat <- do.call(cbind, coef_list)
  
  adj_mat <- .symmetrize(coef_mat)
  adj_mat[which(abs(adj_mat) >= tol)] <- 1
  adj_mat[which(abs(adj_mat) <= tol)] <- 0
  adj_mat
}

.compute_reg_coefficients_cv <- function(dat, lambda = "lambda.1se", verbose = F){
  d <- ncol(dat)

  func <- function(x){
    if(verbose & x %% floor(d/10) == 0) cat('*')
    
    res <- glmnet::cv.glmnet(x = dat[,-x], y = dat[,x], intercept = F)
    vec <- rep(0, d)
    vec[-x] <- as.numeric(stats::coef(res, s = lambda))[-1]
    vec
  }
  
  i <- 0 #debugging purposes only
  foreach::"%dopar%"(foreach::foreach(i = 1:d), func(i))
}

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}

.symmetrize <- function(mat){
  stopifnot(ncol(mat) == nrow(mat))
  (mat + t(mat))/2
}
