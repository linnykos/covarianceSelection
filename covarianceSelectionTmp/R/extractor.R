.split_name <- function(str){
  id <- unlist(strsplit(str,"\\."))[1]
  subregion <- unlist(strsplit(str,"\\."))[2]

  dat <- covarianceSelection::region_subregion
  if(subregion %in% dat$region){
    region <- subregion
  } else {
    stopifnot(subregion %in% dat$subregion)
    region <- dat$region[which(dat$subregion == subregion)]
  }

  dat <- covarianceSelection::brainspan_id
  stopifnot(id %in% dat$Braincode)
  time <- dat$Stage[which(dat$Braincode == id)]

  c(id, region, time)
}

#' Split the dataset into smaller datasets
#'
#' According to the naming convention set in \code{.split_name}
#' and \code{region_subregion}, split the dataset according to the
#' rows The output is a list where each data frame originated
#' from the same ID and same brain region.
#'
#' @param dat data frame
#'
#' @return list of data frames
#' @export
extractor <- function(dat){
  stopifnot(is.data.frame(dat))
  id <- sapply(rownames(dat), .split_name)
  id <- as.factor(apply(id, 2, function(x){paste0(x[1], ".", x[2], ".", x[3])}))

  split(dat, f = id)
}
