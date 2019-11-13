#' Screening data
#'
#' We use a screening method where genes with \code{pv} lower than
#' \code{p_thres} are primary genes, and the secondary genes are selected based on having the having
#' highest
#' (absolute) Pearson correlation with any of the primary genes
#' such that there are \code{num_genes} total.
#'
#' @param dat the nxd matrix
#' @param pv vector of d pvalues
#' @param p_thres numeric
#' @param num_genes numeric
#'
#' @return list of primary and secondary genes, each a vector of integers
#' between 1 and \code{ncol(dat)}
#' @export
screen <- function(dat, pv, p_thres = 0.1, num_genes = 3500){
  stopifnot(ncol(dat) == length(pv))
  stopifnot(is.numeric(dat), is.matrix(dat))

  n <- nrow(dat); d <- ncol(dat)
  
  cor_mat <- abs(stats::cor(dat))

  primary <- which(pv < p_thres)
  # idx <- which(colSums(cor_mat[primary, primary] > cor_thres) > 0)
  # primary <- primary[idx]

  candidates <- c(1:d)[-primary]
  if(length(candidates) == 0 | length(primary) >= num_genes) {
    return(list(primary = primary, secondary = NA, cor_thres = NA))
  }
  
  #find the largest correlation of secondary genes to primary genes
  cor_vec <- apply(cor_mat[candidates, primary], 1, max)
  # secondary <- candidates[which(cor_vec >= cor_thres)]
  idx <- order(cor_vec, decreasing = T)[1:(num_genes - length(primary))]
  cor_thres <- min(cor_vec[idx])
  secondary <- candidates[idx]
  
  stopifnot(length(intersect(primary, secondary)) == 0)

  list(primary = sort(primary), secondary = sort(secondary), cor_thres = cor_thres)
}
