#' Matching
#'
#' Outputs a result so that \code{vec1 == vec2[matching(vec1, vec2)]}.
#'
#' @param vec1 numeric vector
#' @param vec2 numeric vector
#'
#' @return vector of indices
#' @export
matching <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  stopifnot(length(vec1) == length(unique(vec1)))
  stopifnot(length(vec2) == length(unique(vec2)))
  stopifnot(is.vector(vec1), is.vector(vec2))
  stopifnot(all(vec2 %in% vec1), all(vec1 %in% vec2))

  order(vec2)[rank(vec1)]
}
