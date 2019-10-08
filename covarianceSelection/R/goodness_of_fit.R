#' Goodness of fit diagnostic
#'
#' Partitions a list of datasets into two halves, and runs the Cai covariance
#' test on the two halves. This function first centers each half before
#' aggregating. The function tries different unique ways to partition the data
#' into two halves, and returns a p-value for each attempted partition.
#'
#' @param dat_list list of data matrices with same number of columns
#' @param num_partitions number of partitions to attempt
#' @param trials number of trials
#' @param cores number of cores
#' @param verbose boolean
#'
#' @return vector of p-values
#' @export
goodness_of_fit <- function(dat_list, num_partitions, trials = 500, cores = 1,
                            verbose = F){
  dat_list <- lapply(dat_list, function(x){
    scale(x, center = T, scale = F)
  })

  partitions <- .partition_data(length(dat_list), num_partitions)
  sapply(1:length(partitions), function(x){
    if(verbose && length(partitions) > 10 && x %% floor(length(partitions)/10) == 0) cat('*')
    dat1 <- do.call(rbind, dat_list[partitions[[x]][[1]]])
    dat2 <- do.call(rbind, dat_list[partitions[[x]][[2]]])
    cai_test(dat1, dat2, trials = trials, cores = cores)
  })
}

#' Enumerate partitions
#'
#' @param n positive integer
#' @param num_partitions number of partitions
#'
#' @return list of partitions, each a list of 2 vector of integers
.partition_data <- function(n, num_partitions){
  max_val <- 2^(n-1) - 1
  if(num_partitions >= max_val){
    idx <- c(1:max_val)
  } else if(max_val < 1e6) {
    idx <- sample(1:max_val, num_partitions)
  } else {
    idx <- ceiling(stats::runif(num_partitions) * max_val)
  }

  lis <- lapply(idx, function(x){
    vec <- binaryLogic::as.binary(x)
    if(length(vec) <= n-1){
      vec <- c(rep(0, n-1-length(vec)), vec)
    }

    if(all(vec == 0)){
      list(n, c(1:(n-1)))
    } else {
      list(which(vec == 1), c(which(vec == 0), n))
    }
  })
}
