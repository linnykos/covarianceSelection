#' Stepdown to test all pairs of datasets to see if they have equal covariance
#'
#' @param dat_list list of data matrices with same number of columns
#' @param trials number of trials
#' @param alpha alpha level
#' @param denominator boolean
#' @param cores number of cores
#' @param verbose boolean for verbose
#'
#' @return indices for \code{combn(length(dat_list), 2)} that correspond to the
#' hypotheses that passed
#' @export
stepdown <- function(dat_list, trials = 100, alpha = 0.05, denominator = T, cores = 1,
                     verbose = F){
  doMC::registerDoMC(cores = cores)

  dat_list <- lapply(dat_list, scale, center = T, scale = F)
  len <- length(dat_list)
  combn_mat <- utils::combn(len, 2)
  idx_all <- rep(TRUE, ncol(combn_mat))

  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  cov_list <- lapply(dat_list, function(x){n <- nrow(x); (n-1)/n*stats::cov(x)})
  num_list <- lapply(cov_list, function(x){x[diag_idx]})
  if(denominator){
    denom_list <- .compute_all_denom(dat_list, cov_list)
  } else {
    denom_list <- lapply(1:length(dat_list), function(x){1})
  }

  t_vec <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat,
                                  squared = denominator)

  func <- function(i){
    if(verbose && i %% floor(trials/10) == 0) cat('*')
    set.seed(round*10*i)
    noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
    remaining_pairs <- which(idx_all)
    combn_short <- combn_mat[, remaining_pairs, drop = F]
    if(any(is.na(combn_short))) return(Inf)

    remaining_idx <- unique(as.vector(combn_short))
    num_list <- .compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, diag_idx,
                                                  remaining_idx = remaining_idx)

    if(denominator){
      max(abs(.compute_all_test_stat(num_list, denom_list, combn_mat = combn_short)))
    } else {
      .compute_max_accelerated(num_list, combn_mat = combn_short)
    }
  }

  round <- 1
  while(TRUE){
    if(sum(idx_all) == 0) break()

    i <- 0 #debugging purposes
    t_boot <- as.numeric(unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                            func(i))))
    cutoff <- as.numeric(stats::quantile(abs(t_boot), 1-alpha))
    idx <- intersect(which(abs(t_vec) >= cutoff), which(idx_all))

    if(length(idx) == 0) break()

    idx_all[idx] <- FALSE
    round <- round+1
    if(verbose) print(paste0("In stepdown, finished round ", round, " with ",
                             sum(idx_all), " null hypothesis remaining"))
  }

  which(idx_all)
}

#' Compute all of the denominators
#'
#' @param dat_list list of data matrices
#' @param cov_list list of covariance matrices
#'
#' @return a list of denominator vectors
.compute_all_denom <- function(dat_list, cov_list = as.list(rep(NA, length(dat_list)))){
  lapply(1:length(dat_list), function(x){
    .compute_variance(dat_list[[x]], cov_list[[x]])
  })
}

#########

#' Compute the test statistic
#'
#' @param num_list list of numerator vectors
#' @param denom_list list of denominator vectors
#' @param combn_mat matrix of pairs to test for
#' @param squared boolean on whether or not to square the test statistic
#'
#' @return vector of all the bootstrap statistics
.compute_all_test_stat <- function(num_list, denom_list,
                                   combn_mat = utils::combn(length(num_list), 2),
                                   squared = T){
  stopifnot(length(num_list) == length(denom_list))

  sapply(1:ncol(combn_mat), function(x){
    .compute_covStat(num_list[[combn_mat[1,x]]], num_list[[combn_mat[2,x]]],
                     denom_list[[combn_mat[1,x]]], denom_list[[combn_mat[2,x]]],
                     squared = squared)
  })
}

#' Compute all the numerator bootstraps
#'
#' @param dat_list list of data matrices
#' @param noise_list list of noise vectors
#' @param cov_list list of covariance matrices
#' @param idx vector of the indices of the lower triangle
#' @param remaining_idx indices to compute the numerator of
#'
#' @return a list of vectors
.compute_all_numerator_bootstrap <- function(dat_list, noise_list, cov_list, idx,
                                              remaining_idx){
  k <- length(dat_list)

  lis <- vector("list", k)
  lis[remaining_idx] <- lapply(remaining_idx, function(x){.compute_bootSigma(dat_list[[x]], noise_list[[x]],
                                               cov_list[[x]], idx)})
  lis
}

#############

.shorten_combn <- function(combn_mat, num_in){
  if(ncol(combn_mat) == 1){if(all(combn_mat <= num_in)) return(combn_mat) else return(NA)}

  res <- apply(combn_mat, 1, function(x){x <= num_in})
  idx <- which(apply(res, 1, function(x){all(x)}))
  if(length(idx) == 0) return(NA)
  combn_mat[,idx, drop = F]
}
