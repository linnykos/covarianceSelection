#' Nonparanormal transformation
#'
#' @param dat \code{n} by \code{d} matrix
#' @param density_list list of \code{d} \code{density} objects
#' @param mean_vec numeric vector of length \code{d}
#' @param sd_vec positive numeric vector of length \code{d}
#'
#' @return \code{n} by \code{d} matrix
#' @export
nonparanormal_transformation <- function(dat, density_list,
                                         mean_vec, sd_vec){
  stopifnot(all(sapply(density_list, class) == "density"), length(density_list) == ncol(dat))
  stopifnot(length(mean_vec) == length(sd_vec))
  d <- ncol(dat)
  
  for(i in 1:d){
    dat[,i] <- .variable_transformation(dat[,i], density_list[[i]],
                                        mean_vec[i], sd_vec[i])
  }
  
  dat
}

.variable_transformation <- function(vec, den, mean_val = 0, sd_val = 1){
  val <- stats::pnorm(vec, mean = mean_val, sd = sd_val)
  
  # make a quantile vector of the target density
  quantile_vec <- cumsum(den$y)
  quantile_vec <- quantile_vec/max(quantile_vec)
  
  idx <- sapply(val, function(x){
    tmp <- x - quantile_vec
    tmp[tmp < 0] <- NA
    which.min(tmp)
  }) 
  
  difference <- (val - quantile_vec[idx])/(quantile_vec[idx+1]-quantile_vec[idx])
  
  den$x[idx] + difference*(den$x[idx+1]-den$x[idx])
}