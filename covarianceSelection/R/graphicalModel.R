#' Graphical model estimate
#'
#' Estimated using neighbhorhood selection, cross validation to select lambda
#'
#' @param dat the nxd matrix
#' @param verbose boolean
#'
#' @return a precision matrix estimate, dxd
#' @export
graphicalModel <- function(dat, lambda = "lambda.1se", verbose = F){
  n <- nrow(dat); d <- ncol(dat)

  if(verbose) print("Starting to estimate coefficients")
  coef_list <- .compute_reg_coefficients_cv(dat, lambda = lambda, verbose = verbose)
  coef_mat <- do.call(cbind, coef_list)
  
  if(verbose) print("Starting to estimate sigma")
  sigma_vec <- .compute_sigma_vec(dat, coef_mat)

  if(verbose) print("Finalizing")
  prec_mat <- sapply(1:d, function(x){
    vec <- numeric(d)
    vec<- -coef_mat[,x]/sigma_vec[x]
    vec[x] <- 1/sigma_vec[x]

    vec
  })

  .symmetrize(prec_mat)
}

.compute_reg_coefficients_cv <- function(dat, lambda = "lambda.1se", verbose = F){
  d <- ncol(dat)

  func <- function(x){
    if(verbose & x %% floor(d/10) == 0) cat('*')
    
    res <- glmnet::cv.glmnet(x = dat[,-x], y = dat[,x], intercept = F)
    vec <- rep(0, d)
    vec[-x] <- as.numeric(stats::coef(res, s = "lambda.1se"))[-1]
    vec
  }
  
  i <- 0 #debugging purposes only
  foreach::"%dopar%"(foreach::foreach(i = 1:d), func(i))
}

.compute_sigma_vec <- function(dat, coef_mat, cores = 1){
  doMC::registerDoMC(cores = cores)
  d <- ncol(dat); n <- nrow(dat)

  func <- function(x){
    coef_vec <- coef_mat[,x]; coef_vec <- coef_vec[-x]
    .l2norm(dat[,x] - dat[,-x]%*%coef_vec)^2/n
  }

  i <- 0 #debugging purposes only
  as.numeric(unlist(foreach::"%dopar%"(foreach::foreach(i = 1:d), func(i))))
}

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}

.symmetrize <- function(mat, bool_eigen = F){
  stopifnot(ncol(mat) == nrow(mat))

  mat <- (mat + t(mat))/2

  if(bool_eigen){
    res <- eigen(mat)
    res$values[res$values < 0] <- 0

    res$vectors %*% diag(res$values) %*% t(res$vectors)
  } else{
    mat
  }

}
