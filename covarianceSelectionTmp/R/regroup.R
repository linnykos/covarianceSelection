#' Regroup a list of matrices into another list of matrices
#'
#' \code{dat_list} is a list of matrices, and \code{idx_list} is
#' a list of integer vectors, where each integer from 1 to \code{length(dat_list)}
#' appears exactly once. Outputs a list of matrices of length \code{length(idx_list)}
#' matrices.
#'
#' @param dat_list a list of matrices
#' @param idx_list a list of integer vectors
#'
#' @return a list of matrices
#' @export
regroup <- function(dat_list, idx_list){
  .is.listOfMatrix(dat_list, "dat_list")
  .is.listofNumeric(idx_list, "idx_list")
  .regroup_check(length(dat_list), idx_list)

  dat <- do.call(cbind, dat_list)
  idx.vec <- .populate(sapply(dat_list, ncol), idx_list)

  lis <- vector("list", length(idx_list))

  for(i in 1:length(idx_list)){
    idx <- which(idx.vec == i)
    lis[[i]] <- dat[,idx]
  }

  names(lis) <- names(idx_list)
  lis
}

.populate <- function(ncol.vec, idx_list){
  idx.vec <- .partition_renamed(idx_list)
  rep(idx.vec, times = ncol.vec)
}

.partition_renamed <- function(idx_list){
  idx.vec <- unlist(idx_list)
  n <- max(idx.vec)

  vec <- rep(0, n)

  for(i in 1:length(idx_list)){
    vec[idx_list[[i]]] <- i
  }

  vec
}

.regroup_check <- function(len, idx_list){
  vec <- unlist(idx_list)

  .is.nonNegInteger(vec)
  if(length(vec) != len) stop(paste("idx_list is missing elements",
    "compared to dat_list"))
  if(!all(sort(vec) == 1:len)) stop(paste("idx_list does not have",
    "all consecutive integers"))

  TRUE
}
