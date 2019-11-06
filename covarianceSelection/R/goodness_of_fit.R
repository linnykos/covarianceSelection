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
goodness_of_fit <- function(dat, permutations = 250, trials = 100, prob = 1, verbose = F){
  n <- nrow(dat)
  
  sapply(1:permutations, function(x){
    if(verbose && trials > 10 && x %% floor(trials/10) == 0) cat('*')
    split1 <- sample(1:n, round(n/2))
    
    dat1 <- dat[split1,]
    dat2 <- dat[-split1,]
    cai_test(dat1, dat2, trials = trials, prob = prob)
  })
}