#' Goodness of fit diagnostic
#'
#' Partitions a list of datasets into two halves, and runs the Cai covariance
#' test on the two halves. This function first centers each half before
#' aggregating. The function tries different unique ways to partition the data
#' into two halves, and returns a p-value for each attempted partition.
#'
#' @param dat_list list of data matrices with same number of columns
#' @param permutations number of permutations to run
#' @param trials number of trials
#' @param prob \code{prob} parameter for \code{cai_test}
#' @param verbose boolean
#'
#' @return vector of p-values
#' @export
goodness_of_fit <- function(dat_list, permutations = 250, trials = 100, prob = 1, 
                            verbose = F){
  n <- length(dat_list)
  
  sapply(1:permutations, function(x){
    if(verbose && permutations > 10 && x %% floor(permutations/10) == 0) cat('*')
    split1 <- sample(1:n, round(n/2))
    
    dat1 <- do.call(rbind, dat_list[split1])
    dat2 <- do.call(rbind, dat_list[-split1])
    cai_test(dat1, dat2, trials = trials, prob = prob)
  })
}