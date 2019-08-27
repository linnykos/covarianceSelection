roc_region <- function(tpr_mat, fpr_mat,
                       fpr_seq = seq(0, 1, length.out = 21),
                       quant = c(0.1, 0.9)){
  stopifnot(nrow(tpr_mat) == nrow(fpr_mat),
            ncol(tpr_mat) == ncol(fpr_mat))
  stopifnot(all(tpr_mat >= 0), all(tpr_mat <= 1),
            all(fpr_mat >= 0), all(fpr_mat <= 1))
  
  res <- lapply(fpr_seq, function(val){
    index_list <- .find_inbetween_indices(fpr_mat, val)
    percentage_list <- .compute_percentage_inbetween(fpr_mat, val, index_list)
    tpr_val <- .compute_y_coordinate(tpr_mat, index_list, percentage_list)
    
    c(stats::quantile(sapply(tpr_val, min), probs = quant[1]), 
      stats::median(unlist(tpr_val)),
      stats::quantile(sapply(tpr_val, max), probs = quant[2]))
  })
  
  res <- do.call(rbind, res)
  
  list(lower_tpr = res[,1], median_tpr = res[,2], upper_tpr = res[,3])
}

.find_inbetween_indices <- function(fpr_mat, val){
  trials <- ncol(fpr_mat)
  n <- nrow(fpr_mat)
  
  res <- vector("list", trials)
  for(j in 1:trials){
    for(i in 1:(n-1)){
      if(sign(fpr_mat[i,j] - val) != sign(fpr_mat[i+1,j] - val)){
        res[[j]] <- c(res[[j]], i)
      }
    }
  }
  
  res
}

.compute_percentage_inbetween <- function(fpr_mat, val, index_list){
  trials <- ncol(fpr_mat)
  n <- nrow(fpr_mat)
  
  lapply(1:trials, function(j){
    len <- length(index_list[[j]])
    tmp <- sapply(1:len, function(x){
      val1 <- fpr_mat[index_list[[j]][x]+1,j]
      val2 <- fpr_mat[index_list[[j]][x],j]
      stopifnot(sign(val1 - val) != sign(val2 - val))
      
      (val - val1)/(val2 - val1)
    })
  })
}

.compute_y_coordinate <- function(tpr_mat, index_list, percentage_list){
  trials <- ncol(tpr_mat)
  n <- nrow(tpr_mat)
  
  lapply(1:trials, function(j){
    len <- length(index_list[[j]])
    tmp <- sapply(1:len, function(x){
      val1 <- tpr_mat[index_list[[j]][x]+1,j]
      val2 <- tpr_mat[index_list[[j]][x],j]
      
      val1+(val2-val1)*percentage_list[[j]][x]
    })
  })
}