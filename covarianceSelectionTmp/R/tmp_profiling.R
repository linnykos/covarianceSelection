.compute_bootSigma_tmp <- function(mat, noise_vec, cov_mat){
  n <- nrow(mat)
  mat <- scale(mat, center = TRUE, scale = FALSE)
  t(mat)%*%diag(noise_vec/n)%*%mat - (sum(noise_vec)/n)*cov_mat
}

#' Something
#'
#' @param dat_list Something
#' @param noise_list Something
#' @param cov_list Something
#' @param remaining_idx Something
#'
#' @return Something
#' @export
compute_all_numerator_bootstrap_tmp1 <- function(dat_list, noise_list, cov_list,
                                                 remaining_idx){
  k <- length(dat_list)
  
  lis <- vector("list", k)
  lis[remaining_idx] <- lapply(remaining_idx, function(x){.compute_bootSigma_tmp(dat_list[[x]], noise_list[[x]],
                                                                                 cov_list[[x]])})
  lis
}

#' Something
#'
#' @param dat_list Something
#' @param noise_list Something
#' @param cov_list Something
#' @param remaining_idx Something
#'
#' @return Something
#' @export
compute_all_numerator_bootstrap_tmp2 <- function(dat_list, noise_list, cov_list,
                                                 remaining_idx){
  c_compute_all_numerator_bootstrap(dat_list, noise_list, cov_list, remaining_idx)
}