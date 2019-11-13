#' Stepdown path
#'
#' @param dat_list list of data matrices with same number of columns
#' @param trials number of trials
#' @param iterations number of iterations
#' @param cores number of cores
#' @param verbose boolean for verbose
#' @param file temporary save filepath
#' @param prob numeric between 0 and 1
#'
#' @return a \code{stepdown} object
#' @export
stepdown_path <- function(dat_list, trials = 100, iterations = 15, cores = 1,
                          verbose = F, file = NA, prob = 1, squared = T){
  doMC::registerDoMC(cores = cores)
  boot_list <- vector("list", length = iterations)
  if(!is.na(file)) save(boot_list, file = file)
  
  if(verbose)  print(paste0("Entered stepdown function: ", Sys.time()))

  dat_list <- lapply(dat_list, scale, center = T, scale = F)
  diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
  len <- length(dat_list)
  combn_mat <- utils::combn(len, 2)
  
  num_list <- lapply(dat_list, function(x){.compute_sigma(x, diag_idx)})
  if(squared){
    denom_list <- .compute_all_denom(dat_list, num_list, diag_idx)
  } else {
    denom_list <- lapply(1:length(dat_list), function(x){1})
  }

  t_vec <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat, squared = squared,
                                  prob = prob)

  if(verbose)  print(paste0("Starting to run heavy parallel computation: ", Sys.time()))
  
  func <- function(i, round){
    if(verbose && i %% floor(trials/10) == 0) cat('*')
    set.seed((round-1)*trials + i)
    noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
    boot_num_list <- .compute_all_numerator_bootstrap(dat_list, noise_list, num_list, diag_idx,
                                                 remaining_idx = 1:len)

    .compute_all_test_stat(boot_num_list, denom_list, combn_mat = combn_mat, squared = squared, prob = 1)
  }

  for(x in 1:iterations){
    if(verbose) print(paste0("On iteration ", x))
    i <- 0 #debugging purposes
    
    boot_list[[x]] <- do.call(rbind, foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i, round = x)))
    res <- list(t_vec = t_vec, boot = boot_list)
    if(!is.na(file)) save(res, file = file)
  }

  structure(list(t_vec = t_vec, boot = boot_list), class = "stepdown")
}

#' Stepdown choose
#'
#' @param stepdown_obj stepdown object
#' @param alpha alpha level
#' @param return_pvalue boolean if the naive p-values should be returned
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
    
    cutoff <- stats::quantile(abs(t_boot), 1-alpha)
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
      # aggregate all the potential test statistics
      test_boot_vec <- unlist(lapply(1:length(stepdown_obj$boot), function(k){
        stepdown_obj$boot[[k]][,i]
      }))
      
      length(which(abs(test_boot_vec) > abs(stepdown_obj$t_vec[i])))/length(test_boot_vec)
    })
  } else {
    pval <- NA
  }

  list(null_idx = which(idx_all), pval = pval)
}
