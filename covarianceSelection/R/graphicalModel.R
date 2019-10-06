#' Graphical model estimate
#'
#' Estimated using neighbhorhood selection, cross validation to select lambda
#'
#' @param dat the nxd matrix
#' @param verbose boolean
#'
#' @return a precision matrix estimate, dxd
#' @export
graphicalModel <- function(dat, primary_idx, lambda = "lambda.1se", verbose = F, tol = 1e-6){
  n <- nrow(dat); d <- ncol(dat)

  if(verbose) print("Starting to estimate coefficients")
  coef_list <- .compute_reg_coefficients_cv(dat, primary_idx, lambda = lambda, verbose = verbose)
  coef_mat <- sapply(coef_list, function(x){ x$vec })
  lambda_vec <- sapply(coef_list, function(x){ x$lambda })
  if(abs(diff(range(lambda_vec))) <= tol) lambda_vec <- min(lambda_vec)
  stopifnot(nrow(coef_mat) == ncol(dat), ncol(coef_mat) == length(primary_idx))
  
  adj_mat <- cbind(coef_mat, matrix(0, nrow = nrow(coef_mat), ncol = ncol(coef_mat) - nrow(coef_mat)))
  adj_mat <- .symmetrize(coef_mat)
  adj_mat[which(abs(adj_mat) >= tol)] <- 1
  adj_mat[which(abs(adj_mat) <= tol)] <- 0
  adj_mat <- Matrix::Matrix(adj_mat, sparse = T)
  
  list(adj_mat = adj_mat, lambda_vec = lambda_vec)
}

graphicalModel_range <- function(dat,  primary_idx, lambda_min, lambda_max, lambda_length = 15, verbose = F, tol = 1e-6){
  lambda_seq <- exp(seq(log(lambda_min), log(lambda_max), length.out = lambda_length))
  
  lapply(lambda_seq, function(x){
    if(verbose) print(x)
    graphicalModel(dat,  primary_idx, lambda = x, verbose = verbose, tol = tol)
  })
}

##################

.compute_reg_coefficients_cv <- function(dat, primary_idx = 1:ncol(dat), lambda = "lambda.1se", verbose = F){
  d <- ncol(dat)

  func <- function(i){
    if(verbose & i %% floor(length(primary_idx)/10) == 0) cat('*')
    
    x <- primary_idx[i]
    if(is.numeric(lambda)){
      res <- glmnet::glmnet(x = dat[,-x], y = dat[,x], intercept = F, lambda = lambda)
      vec <- rep(0, d)
      vec[-x] <- as.numeric(res$beta)
      
      list(vec = vec, lambda = lambda)
    } else {
      res <- glmnet::cv.glmnet(x = dat[,-x], y = dat[,x], intercept = F)
      vec <- rep(0, d)
      vec[-x] <- as.numeric(stats::coef(res, s = lambda))[-1]
      
      list(vec = vec, lambda = res[[which(names(res) == lambda)]])
    }
    
  }
  
  i <- 0 #debugging purposes only
  foreach::"%dopar%"(foreach::foreach(i = 1:length(primary_idx)), func(i))
}

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}

.symmetrize <- function(mat){
  stopifnot(ncol(mat) == nrow(mat))
  (mat + t(mat))/2
}
