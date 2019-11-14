#' Stepdown to test all pairs of datasets to see if they have equal covariance
#'
#' @param dat_list list of data matrices with same number of columns
#' @param trials number of trials
#' @param alpha alpha level
#' @param return_pvalue boolean if the naive p-values should be returned
#' @param cores number of cores
#' @param prob numeric between 0 and 1
#' @param verbose boolean for verbose
#' @param only_test_stat boolean if only the test statistic should be return
#'
#' @return if \code{only_test_stat} is \code{FALSE}:
#' a list containing \code{null_idx} for indices for \code{combn(length(dat_list), 2)} that correspond to the
#' hypotheses that passed and \code{pval} (which will be \code{NA} if \code{return_pvalue} is
#' \code{FALSE}) or the p-values
#' @export
stepdown <- function(dat_list, trials = 100, alpha = 0.05, return_pvalue = F, cores = 1,
                     prob = 1, verbose = F, only_test_stat = F, squared = T){
  doMC::registerDoMC(cores = cores)
  
  if(verbose)  print(paste0("Entered stepdown function: ", Sys.time()))

  dat_list <- lapply(dat_list, scale, center = T, scale = F)
  len <- length(dat_list)
  combn_mat <- utils::combn(len, 2)
  idx_all <- rep(TRUE, ncol(combn_mat))

  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  num_list <- lapply(dat_list, function(x){.compute_sigma(x, diag_idx)})
  denom_list <- .compute_all_denom(dat_list, num_list, diag_idx)

  t_vec <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat, squared = squared, prob = prob)
  if(only_test_stat) return(t_vec)
  
  if(verbose)  print(paste0("Starting to run heavy parallel computation: ", Sys.time()))

  func <- function(i){
    if(verbose && i %% floor(trials/10) == 0) cat('*')
    set.seed((round-1)*trials + i)
    noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
    
    remaining_pairs <- which(idx_all)
    combn_short <- combn_mat[, remaining_pairs, drop = F]
    if(any(is.na(combn_short))) return(Inf)

    remaining_idx <- unique(as.vector(combn_short))
    boot_num_list <- .compute_all_numerator_bootstrap(dat_list, noise_list, num_list, diag_idx,
                                                  remaining_idx = remaining_idx)

    boot_t_vec <- .compute_all_test_stat(boot_num_list, denom_list, combn_mat = combn_short,
                                         prob = 1, squared = squared)
    
    if(round == 1 & return_pvalue){
      list(val = max(abs(boot_t_vec)), boot_t_vec = boot_t_vec)
    } else {
      list(val = max(abs(boot_t_vec)), boot_t_vec = NA)
    }
  }

  round <- 1
  round_1_boot_t_vec <- NA
  
  while(TRUE){
    if(sum(idx_all) == 0) break()

    i <- 0 #debugging purposes
    res <- foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i))
    t_boot <- sapply(res, function(x){x$val})
    if(round == 1 & return_pvalue) round_1_boot_t_vec <- sapply(res, function(x){x$boot_t_vec})
    
    cutoff <- stats::quantile(abs(t_boot), 1-alpha)
    idx <- intersect(which(abs(t_vec) >= cutoff), which(idx_all))

    if(length(idx) == 0) break()

    idx_all[idx] <- FALSE
    round <- round + 1
    if(verbose) print(paste0("In stepdown, finished round ", round, " with ",
                             sum(idx_all), " null hypothesis remaining"))
  }

  if(return_pvalue){
    stopifnot(length(t_vec) == nrow(round_1_boot_t_vec))
    pval <- sapply(1:length(t_vec), function(i){
      length(which(abs(round_1_boot_t_vec[i,]) > abs(t_vec[i])))/ncol(round_1_boot_t_vec)
    })
  } else {
    pval <- NA
  }
  
  list(null_idx = which(idx_all), pval = pval)
}

#' Compute all of the denominators
#'
#' @param dat_list list of data matrices
#' @param num_list list of covariance matrices
#' @param idx vector of the indices of the lower triangle
#'
#' @return a list of denominator vectors
.compute_all_denom <- function(dat_list, num_list, idx){
  lapply(1:length(dat_list), function(x){
    .compute_variance(dat_list[[x]], num_list[[x]], idx)
  })
}

#########

#' Compute the test statistic
#'
#' @param num_list list of numerator vectors
#' @param denom_list list of denominator vectors
#' @param combn_mat matrix of pairs to test for
#' @param squared boolean on whether or not to square the test statistic
#' @param prob numeric between 0 and 1
#'
#' @return vector of all the bootstrap statistics
.compute_all_test_stat <- function(num_list, denom_list,
                                   combn_mat = utils::combn(length(num_list), 2),
                                   squared = T, prob = 1){
  stopifnot(length(num_list) == length(denom_list))

  sapply(1:ncol(combn_mat), function(x){
    .compute_covStat(num_list[[combn_mat[1,x]]], num_list[[combn_mat[2,x]]],
                     denom_list[[combn_mat[1,x]]], denom_list[[combn_mat[2,x]]],
                     squared = squared, prob = prob)
  })
}

#' Compute all the numerator bootstraps
#'
#' @param dat_list list of data matrices
#' @param noise_list list of noise vectors
#' @param num_list list of covariance matrices
#' @param idx vector of the indices of the lower triangle
#' @param remaining_idx indices to compute the numerator of
#'
#' @return a list of vectors
.compute_all_numerator_bootstrap <- function(dat_list, noise_list, num_list, idx,
                                              remaining_idx){
  k <- length(dat_list)

  lis <- vector("list", k)
  lis[remaining_idx] <- lapply(remaining_idx, function(x){.compute_bootSigma(dat_list[[x]], noise_list[[x]],
                                                                             num_list[[x]], idx)})
  lis
}