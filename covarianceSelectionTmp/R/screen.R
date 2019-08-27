#' Screening data
#'
#' We use a screening method where genes with \code{pv} lower than
#' \code{pthres} are primary genes. \code{num_genes} are selected, where
#' the remaining secondary genes are selected based on having the having
#' (absolute) Pearson correlation with any of the primary genes.
#'
#' @param dat the nxd matrix
#' @param pv vector of d pvalues
#' @param pthres numeric
#' @param num_genes numeric
#'
#' @return list of primary and secondary genes, each a vector of integers
#' between 1 and \code{ncol(dat)}
#' @export
screen <- function(dat, pv, pthres = 0.1, num_genes = min(100, ceiling(nrow(dat)/10))){
  stopifnot(ncol(dat) == length(pv))
  stopifnot(is.numeric(dat), is.matrix(dat))

  n <- nrow(dat); d <- ncol(dat)

  primary <- which(pv < pthres)

  cor_mat <- abs(stats::cor(dat))
  candidates <- c(1:d)[-primary]
  if(length(candidates) == 0 | num_genes <= length(primary)) {
    return(list(primary = primary, secondary = NA))
  }

  #find the largest correlation of secondary genes to primary genes
  cor_vec <- apply(cor_mat[candidates, primary], 1, max)
  secondary <- candidates[order(cor_vec, decreasing = T)[1:(num_genes - length(primary))]]

  list(primary = sort(primary), secondary = sort(secondary))
}
