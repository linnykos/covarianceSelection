#' Binning of subjects
#'
#' This must be the in the format of ID.REGION.(whatever), for example,
#' "HSB100.PFC.6", where the ID is prefixed with "HSB". The splits are
#' a vector of integers between 1 and 15 in increasing order.
#'
#' @param vec a vector of valid subjects
#' @param splits a vector of integers
#'
#' @return a table
#' @export
binning <- function(vec, splits = c(2,5,8,15)){
  stopifnot(all(diff(splits) > 0), min(splits) >= 1, max(splits) <= 15)

  name_mat <- sapply(vec, .split_name)
  time <- .time_to_splits(as.numeric(name_mat[3,]), splits = splits)

  region <- factor(as.character(name_mat[2,]), levels = c("PFC", "VIIAS", "SHA", "MDCBC"))

  table(region, time)
}

.time_to_splits <- function(vec, splits = c(2,5,8,15)){
  vec[vec <= splits[1]] <- 1
  for(j in 2:length(splits)){
    vec[intersect(which(vec >= splits[j-1]+1), which(vec <= splits[j]))] <- j
  }

  factor(vec, levels = c(1,2,3,4))
}
