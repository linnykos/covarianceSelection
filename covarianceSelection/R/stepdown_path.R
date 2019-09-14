#' Stepdown path
#'
#' @param dat_list list of data matrices with same number of columns
#' @param trials number of trials
#' @param iterations number of iterations
#' @param cores number of cores
#' @param verbose boolean for verbose
#'
#' @return a \code{stepdown} object
#' @export
stepdown_path <- function(dat_list, trials = 100, iterations = 15, cores = 1,
                          verbose = F){
  doMC::registerDoMC(cores = cores)

  dat_list <- lapply(dat_list, scale, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  len <- length(dat_list)
  combn_mat <- utils::combn(len, 2)
  
  num_list <- lapply(dat_list, function(x){.compute_sigma(x, diag_idx)})
  denom_list <- .compute_all_denom(dat_list, num_list, diag_idx)

  t_vec <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat, squared = T)

  if(verbose)  print(paste0("Starting to run heavy parallel computation: ", Sys.time()))
  
  func <- function(i, round){
    if(verbose && i %% floor(trials/10) == 0) cat('*')
    set.seed(round*10*i)
    noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
    boot_num_list <- .compute_all_numerator_bootstrap(dat_list, noise_list, num_list, diag_idx,
                                                 remaining_idx = 1:len)

    .compute_all_test_stat(boot_num_list, denom_list, combn_mat = combn_mat)
  }

  boot_list <- lapply(1:iterations, function(x){
    if(verbose) print(paste0("On iteration ", x))
    i <- 0 #debugging purposes
    
    do.call(rbind, foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i, round = x)))
  })

  structure(list(t_vec = t_vec, boot = boot_list), class = "stepdown")
}

#' Stepdown choose
#'
#' @param stepdown_obj stepdown object
#' @param alpha alpha level
#' @param verbose boolean for verbose
#'
#' @return indices for \code{combn(length(dat_list), 2)} that correspond to the
#' hypotheses that passed
#' @export
stepdown_choose <- function(stepdown_obj, alpha = 0.05, return_pvalue = F, verbose = F){
  stopifnot(class(stepdown_obj) == "stepdown")
  len <- length(stepdown_obj$t_vec)
  idx_all <- rep(TRUE, len)

  round <- 1
  while(TRUE){
    if(sum(idx_all) == 0) break()

    boot_mat <- stepdown_obj$boot[[round]][, which(idx_all), drop = F]
    t_boot <- apply(boot_mat, 1, function(x){abs(max(x))})
    
    cutoff <- stats::quantile(t_boot, 1-alpha)
    idx <- intersect(which(abs(stepdown_obj$t_vec) >= cutoff), which(idx_all))

    if(length(idx) == 0) break()

    idx_all[idx] <- FALSE
    if(verbose) print(paste0("In stepdown, finished round ", round, " with ",
                             sum(idx_all), " null hypothesis remaining"))
    round <- round+1

    if(round > length(stepdown_obj$boot)) {
      warning("Function ran out of bootstrap rounds to work with.")
      break()
    }
  }
  
  if(return_pvalue){
    pval <- sapply(1:length(stepdown_obj$t_vec), function(i){
      length(which(abs(stepdown_obj$boot[[1]][,i]) > abs(stepdown_obj$t_vec[i])))/nrow(stepdown_obj$boot[[1]])
    })
  } else {
    pval <- NA
  }

  list(null_idx = which(idx_all), pval = pval)
}
