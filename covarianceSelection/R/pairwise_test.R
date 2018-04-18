#' Pairwise Test Function
#'
#' Uses Cai's covariance test between all the elements in a list.
#' Only the lower-triangular values are filled in.
#'
#' If \code{dat_list2}
#' is provided, then instead this function apply's Cai's test
#' between every element in \code{dat_list} and \code{dat_list2}. In
#' this case, the returned matrix has \code{dat_list} elements
#' along the columns, and \code{dat_list2} elements along the rows.
#'
#' @param dat_list list of matrices
#' @param dat_list2 (optional) list of matrices
#' @param trials number of trials
#' @param quant quantile
#' @param cores number of cores
#'
#' @return matrix of p-values
#' @export
pairwise_test <- function(dat_list, dat_list2 = NA, trials = 100, quant = 1,
                                 cores = 1){
  stopifnot(is.list(dat_list))
  stopifnot(any(sapply(dat_list, is.matrix)), any(sapply(dat_list, is.numeric)),
            length(unique(sapply(dat_list, ncol))) == 1)
  if(!any(is.na(dat_list2))){
    stopifnot(is.list(dat_list2))
    stopifnot(any(sapply(dat_list2, is.matrix)), any(sapply(dat_list2, is.numeric)),
              length(unique(sapply(dat_list2, ncol))) == 1)
    stopifnot(ncol(dat_list) == ncol(dat_list2))
  }

  doMC::registerDoMC(cores = cores)

  n <- length(dat_list)
  i <- 1 #debugging purposes

  if(any(is.na(dat_list2))){
    comb_mat <- utils::combn(n, 2)
    p_vec <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(comb_mat)),
                                       cai_test(dat_list[[comb_mat[1,i]]],
                                                dat_list[[comb_mat[2,i]]],
                                                trials, quant)))

    p_mat <- .convert_combn2mat(p_vec, n)
    rownames(p_mat) <- names(dat_list)
  } else {
    m <- length(dat_list2)
    comb_mat <- rbind(rep(1:n, each = m), rep(1:m, times = n))
    p_vec <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(comb_mat)),
                                       cai_test(dat_list[[comb_mat[1,i]]],
                                                dat_list2[[comb_mat[2,i]]],
                                                trials, quant)))
    p_mat <- matrix(p_vec, nrow = m, ncol = n)
    rownames(p_mat) <- names(dat_list2)
  }

  colnames(p_mat) <- names(dat_list)

  p_mat
}

.convert_combn2mat <- function(vec, d){
  mat <- matrix(0, d, d)
  mat[lower.tri(mat)] <- vec
  mat
}
