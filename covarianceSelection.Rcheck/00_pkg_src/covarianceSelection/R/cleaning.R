#' Gene synonyms
#' 
#' This function is an in-house way to map genes to their most commonly-known form
#' based on the \code{org.Hs.eg.db} package. If there are conflicts (i.e., multiple
#' genes that claim to be the official alias of a gene), we take the most "common" gene (i.e.,
#' how often it appears as the official gene). If there are still ties, we take the 
#' first gene alphabetically.
#' 
#' If no genes are found for a particular gene character, return \code{character(0)}.
#' 
#' At the time of this function's writing, we are using version 3.8.2 of the \code{org.Hs.eg.db} package.
#' 
#' @param vec a character vector
#' @param verbose a boolean
#' 
#' @return a character vec, in order with respect to the input \code{vec}
#' @export
#' 
#' @source The code for this comes from the user Duff at \url{https://www.biostars.org/p/14971/}.
symbol_synonyms <- function(vec, verbose = T){
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
  
  sapply(1:length(vec), function(i){
    if(verbose & i %% max(floor(length(vec)/10),1) == 0) cat('*')
    
    res <- aliasSymbol[which(aliasSymbol[,2] %in% vec[i]), 5]
    if(length(res) == 0) return(NA)
    
    #if there are more than one, take the most common one
    if(length(res) > 1){
      len_vec <- sapply(res, function(x){
        length(which(aliasSymbol[,2] %in% x))
      })
      
      #if there is still a tie, take the first alphabetic one
      res_final <- res[which(len_vec == max(len_vec))]
      if(length(res_final) > 1){
        res_final <- sort(res_final, decreasing = F)[1]
      }
      
      res <- res_final
    }
    
    res
  })
}

average_same_columns <- function(dat){
  col_vec <- colnames(dat)
  tab_vec <- table(col_vec)
  idx <- which(tab_vec > 1)
  remove_idx <- numeric(0)
  
  if(length(idx) > 0){
    for(i in idx){
      col_idx <- which(col_vec == names(tab_vec)[i])
      dat[,col_idx[1]] <- rowMeans(dat[,col_idx])
      remove_idx <- c(remove_idx, col_idx[-1])
    }
    
    dat <- dat[,-remove_idx]
  }
  
  dat
}

####################

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

####################

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
  
  base::split(dat, f = id, drop = F)
}

############################

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
