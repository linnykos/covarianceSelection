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

##########################

#' Report results
#'
#' @param gene vector of characters
#' @param posterior vector of numerics
#' @param pvalue vector of numerics
#' @param Iupdate vector of booleans
#'
#' @return data frame with 4 columns
#' @export
report_results <- function(gene, posterior, pvalue, Iupdate){
  stopifnot(length(gene) == length(posterior), length(gene) == length(pvalue),
            length(gene) == length(Iupdate))
  
  d <- length(posterior)
  rankpost <- sort(posterior)
  localfdr <- sapply(1:d, function(x){mean(rankpost[1:x])})
  
  flocalfdr <- rep(0,d)
  rankp <- rank(posterior, ties.method="random")
  flocalfdr <- localfdr[rankp] #undo the sorting
  
  res = data.frame(gene, pvalue, flocalfdr, Iupdate)
  names(res) = c("Gene","p.value","FDR", "indicator")
  
  res
}

#' Fisher test
#'
#' @param vec1 vector of gene names
#' @param vec2 vector of gene names
#' @param all_vec vector of all gene names
#'
#' @return list containing the contigency matrix and the Fisher p value
#' @export
enrichment_test <- function(vec1, vec2, all_vec){
  mat <- matrix(NA,2,2)
  mat[1,1] <- sum(vec1 %in% vec2)
  mat[1,2] <- length(vec1)-mat[1,1]
  mat[2,1] <- length(vec2)-mat[1,1]
  mat[2,2] <- length(all_vec)-length(which(all_vec %in% c(vec1,vec2)))
  
  list(contigency = mat, fisher = stats::fisher.test(mat))
}


