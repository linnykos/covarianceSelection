#' Graphical model estimate
#'
#' Store the glmnet objects (one for each dimension) in a folder location.
#' Will NOT create the folder, please make sure done beforehand.
#' This function requires cores to be set up for parallelization.
#' The objects are saved as name "res".
#'
#' @param dat the nxd matrix
#' @param filename_func a function to create the file name. it takes in
#' \code{i}, the dimension index, as a paramter
#' @param num_primary integer from 1 to \code{ncol(dat)} that denotes
#' the number of genes to consider
#' @param cores number of cores
#' @param cv boolean on whether or not to cross validate
#'
#' @return void
#' @export
graphicalModel_store <- function(dat, filename_func, num_primary = ncol(dat), cores = 1,
                                 cv = F){
  stopifnot(num_primary >= 1, num_primary <= ncol(dat), num_primary %% 1 == 0)
  n <- nrow(dat); d <- ncol(dat)

  doMC::registerDoMC(cores = cores)
  i <- 1 #debugging purposes
  tmp <- foreach::"%dopar%"(foreach::foreach(i = 1:num_primary),
          .compute_reg_coefficients_save(dat, i, filename_func, cv = cv))

  invisible()
}


.compute_reg_coefficients_save <- function(dat, i, filename_func, cv = F){
  if(cv){
    res <- glmnet::cv.glmnet(x = dat[,-i], y = dat[,i], intercept = F)
    res <- .clean_cv_glmnet(res)
  } else {
    res <- glmnet::glmnet(x = dat[,-i], y = dat[,i], intercept = F)
    res <- .clean_glmnet(res)
  }

  filename <- filename_func(i)
  save(res, file = filename)

  invisible()
}

.clean_glmnet <- function(obj){
  stopifnot(inherits(obj,"glmnet"))

  obj$df <- NULL; obj$dim <- NULL; obj$dev.ratio <- NULL; obj$nulldev <- NULL
  obj$npasses <- NULL; obj$jerr <- NULL; obj$call <- NULL; obj$nobs <- NULL

  rownames(obj$beta) <- numeric(nrow(obj$beta))

  obj
}

.clean_cv_glmnet <- function(obj){
  stopifnot(inherits(obj,"cv.glmnet"))

  obj$lambda <- NULL; obj$cvm <- NULL; obj$cvsd <- NULL; obj$cvup <- NULL
  obj$cvlo <- NULL; obj$nzero <- NULL; obj$name <- NULL

  obj$glmnet.fit <- .clean_glmnet(obj$glmnet.fit)

  obj
}
