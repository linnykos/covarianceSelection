roc_region <- function(tpr_mat, fpr_mat,
                       fpr_seq = seq(0, 1, length.out = 21)){
  stopifnot(nrow(tpr_mat) == nrow(fpr_mat),
            ncol(tpr_mat) == ncol(fpr_mat))
  stopifnot(all(tpr_mat >= 0), all(tpr_mat <= 1),
            all(fpr_mat >= 0), all(fpr_mat <= 1))
  
  
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