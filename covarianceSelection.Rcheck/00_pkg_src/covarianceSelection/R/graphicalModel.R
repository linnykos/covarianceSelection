#' Graphical model estimate
#'
#' Estimated using neighbhorhood selection, cross validation to select lambda
#'
#' @param dat the matrix with \code{n} rows and \code{d} columns
#' @param primary_idx index vector that is a subset of \code{1:ncol(dat)}
#' @param lambda either a character vector (\code{"lambda.1se"} or \code{"lambda.min"}) or a numeric positive scalar
#' @param verbose boolean
#' @param tol numeric
#'
#' @return a list that contains an \code{d} by \code{d} \code{sparseMatrix} encoding the
#' estimated adjacency matrix and a numeric vector \code{lambda_vec}
#' @export
graphicalModel <- function(dat, primary_idx, lambda = "lambda.1se", verbose = F, tol = 1e-6){
  n <- nrow(dat); d <- ncol(dat)

  if(verbose) print("Starting to estimate coefficients")
  coef_list <- .compute_reg_coefficients_cv(dat, primary_idx, lambda = lambda, verbose = verbose)
  coef_mat <- sapply(coef_list, function(x){ x$vec })
  lambda_vec <- sapply(coef_list, function(x){ x$lambda })
  if(abs(diff(range(lambda_vec))) <= tol) lambda_vec <- min(lambda_vec)
  stopifnot(nrow(coef_mat) == ncol(dat), ncol(coef_mat) == length(primary_idx))
  
  adj_mat <- cbind(coef_mat, matrix(0, nrow = nrow(coef_mat), ncol = nrow(coef_mat) - ncol(coef_mat)))
  adj_mat <- .symmetrize(adj_mat)
  adj_mat[which(abs(adj_mat) >= tol)] <- 1
  adj_mat[which(abs(adj_mat) <= tol)] <- 0
  adj_mat <- Matrix::Matrix(adj_mat, sparse = T)
  
  list(adj_mat = adj_mat, lambda_vec = lambda_vec)
}

#' Graphical model estimate for a range of lambda values
#'
#' @param dat the matrix with \code{n} rows and \code{d} columns
#' @param primary_idx index vector that is a subset of \code{1:ncol(dat)}
#' @param lambda_min minimum value of \code{lambda} when using \code{graphicalModel}
#' @param lambda_max maximum value of \code{lambda} when using \code{graphicalModel}
#' @param lambda_length number of \code{lambda} values to try (exponential growth)
#' @param verbose boolean
#' @param tol numeric
#'
#' @return a list, each being an output for \code{graphicalModel} for a different value of \code{lambda}
#' @export
graphicalModel_range <- function(dat,  primary_idx, lambda_min, lambda_max, lambda_length = 15, verbose = F, tol = 1e-6){
  lambda_seq <- seq(lambda_min, lambda_max, length.out = lambda_length)
  
  lapply(lambda_seq, function(x){
    if(verbose) print(x)
    graphicalModel(dat,  primary_idx, lambda = x, verbose = verbose, tol = tol)
  })
}

##################

.compute_reg_coefficients_cv <- function(dat, primary_idx = 1:ncol(dat), lambda = "lambda.1se", verbose = F){
  d <- ncol(dat)
  stopifnot(all(primary_idx %in% 1:d))

  func <- function(i){
    if(verbose & i %% floor(length(primary_idx)/10) == 0) cat('*')
  
    vec <- rep(0, d)
    y_vec <- dat[,i]
    
    if(i %in% primary_idx){
      x_mat <- dat[,-i]; idx_vec <- c(1:d)[-i]
    } else{
      x_mat <- dat[,primary_idx]; idx_vec <- primary_idx
    }
    
    if(is.numeric(lambda)){
      res <- glmnet::glmnet(x = x_mat, y = y_vec, intercept = F, lambda = lambda)
      vec[idx_vec] <- as.numeric(res$beta)
    } else {
      res <- glmnet::cv.glmnet(x = x_mat, y = y_vec, intercept = F)
      vec[idx_vec] <- as.numeric(stats::coef(res, s = lambda))[-1]
    }
    
    list(vec = vec, lambda = ifelse(is.numeric(lambda), lambda, res[[which(names(res) == lambda)]]))
  }
  
  i <- 0 #debugging purposes only
  foreach::"%dopar%"(foreach::foreach(i = primary_idx), func(i))
}

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}

.symmetrize <- function(mat){
  stopifnot(ncol(mat) == nrow(mat))
  (mat + t(mat))/2
}
