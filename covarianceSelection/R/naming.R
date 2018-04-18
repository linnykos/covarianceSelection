#' Gene synonyms
#'
#' @param vec vector of gene names (characters)
#'
#' @return vector of characters
#' @export
symbol_synonyms <- function(vec){
  synonyms <- covarianceSelection::synonyms

  res <- sapply(vec, function(x){synonyms$hash[[x]]})
  res <- sapply(res, function(x){ifelse(is.null(x), NA, x)})

  if(all(is.na(res))) return(rep(NA, length(vec)))
  names(synonyms$syn.list)[res]
}