#' Create filename function
#'
#' If \code{folder} does not exist, create \code{folder}. Returns
#' a function to make the filename. \code{folder} must be a character
#' that ends with "/".
#'
#' @param folder folder name (character)
#' @param prefix prefix that all filenames will receive (character)
#'
#' @return function
#' @export
filename_closure <- function(folder, prefix = ""){
  stopifnot(is.character(folder))
  stopifnot(substr(folder, nchar(folder), nchar(folder)) == "/")

  if(!dir.exists(folder)){
    dir.create(folder, recursive = T)
  }

  func <- function(i){
    paste0(folder, prefix, i, ".RData")
  }
}
