#' Find largest soft clique in graph
#'
#' For a given graph \code{g} (igraph object) and a threshold \code{threshold}
#' (number between 0 and 1), find the largest set of nodes such that
#' the number of edges among this set of nodes in \code{g} is at least
#' \code{threshold} percentage of a full clique.
#'
#' This is returned as a list of indices, so if there are multiple largest sets
#' of the same size, the list will have multiple set of indices.
#'
#' \code{mode} is either "all" or "or". "all" is an approximate form that enables
#' faster computing, but "or" gives exact answers.
#'
#' If \code{target_idx} is not \code{NA}, then maximal cliques are pruned based on
#' their overlap with \code{target_idx} and its length. And after the algorithm,
#' indices in \code{target_idx} are seen if they can be filled in heuristically.
#'
#' @param g \code{igraph} object
#' @param threshold numeric
#' @param mode string
#' @param target_idx vector
#' @param prune boolean
#' @param verbose boolean
#' @param time_limit large number
#' @param max_length maximum length of queue before termination
#'
#' @return list of indices
#' @export
clique_selection <- function(g, threshold = 0.95, mode = "all",
                             target_idx = NA, prune = T,
                             verbose = F, time_limit = 3600,
                             max_length = 50000){
  stopifnot(mode %in% c("all", "or"), length(mode) == 1, class(g) == "igraph")
  
  d <- igraph::vcount(g)
  adj <- as.matrix(igraph::as_adjacency_matrix(g))
  
  clique_list <- lapply(igraph::maximal.cliques(g), function(x){sort(as.numeric(x))})
  if(any(!is.na(target_idx))) clique_list <- .prune_clique(adj, clique_list, target_idx, threshold)
  if(length(clique_list) == 1) return(list(1:ncol(adj)))
  largest_clique <- list(sort(clique_list[[which.max(sapply(clique_list, length))]]))
  
  len <- length(clique_list)
  queue <- .initialize_queue(len)
  hash_history <- hash::hash()
  hash_children <- .initialize_children(clique_list)
  hash_unique <- .initialize_unique(clique_list, d)
  
  start_time <- proc.time()["elapsed"]
  prev_time <- start_time
  stat_vec <- rep(0, 5)
  names(stat_vec) <- c("Length", "Num.checked", "Num.total", "Num.accessed", "Num.max")
  
  while(length(queue) > 0){
    stat_vec["Length"] <- length(queue)
    stat_vec["Num.max"] <- length(largest_clique[[1]])
    
    stopifnot(length(hash_unique) == length(hash_children))
    if(verbose & proc.time()["elapsed"] > prev_time + 5) {
      print(stat_vec)
      prev_time <- proc.time()["elapsed"]
    }
    
    #pop
    string <- dequer::pop(queue)
    pair <- .convert_string_to_pair(string)
    
    #check unique
    node1 <- hash_children[[as.character(pair[1])]]
    node2 <- hash_children[[as.character(pair[2])]]
    element1 <- node1$elements; element2 <- node2$elements
    bool <- .check_nonunique(hash_unique, element1, element2, d)
    if(bool) {hash_history[[string]] <- TRUE; next()}
    
    # check adj
    stopifnot(any(is.null(hash_history[[string]])))
    idx <- sort(unique(c(element1, element2)))
    bool <- .pass_threshold(adj[idx, idx, drop = F], threshold)
    hash_history[[string]] <- bool
    stat_vec["Num.checked"] <- stat_vec["Num.checked"] + 1
    
    # if yes again, add new pairs to queue and to clique_list
    if(bool){
      len <- len+1
      
      if(length(idx) > length(largest_clique[[1]])){
        largest_clique <- list(idx)
      } else if (length(idx) == length(largest_clique[[1]])){
        bool_vec <- sapply(largest_clique, function(x){
          all(idx %in% x)
        })
        
        if(all(!bool_vec)) largest_clique[[length(largest_clique)+1]] <- idx
      }
      
      # update
      stopifnot(length(pair) == 2)
      hash_children[[as.character(len)]] <- .node(elements = idx, children = pair)
      hash_unique[[.convert_indices_to_binary(idx, d)]] <- TRUE
      
      tmp <- .add_to_queue(queue, len, idx, pair, hash_children, hash_unique, hash_history, d)
      stat_vec[c("Num.total", "Num.accessed")] <- stat_vec[c("Num.total", "Num.accessed")] + tmp
    }
    
    if(proc.time()["elapsed"] - start_time > time_limit) {
      warning("Function ran out of time.")
      break()
    }
    
    if(!is.na(max_length) & stat_vec["Length"] > max_length){
      warning("Queue size exceeded limit.")
      break()
    }
  }
  
  if(verbose) print("Exiting function")
  
  if(prune) largest_clique <- .prune_nodes(largest_clique, adj, target_idx)
  
  largest_clique
}

#' Select clique based on intersection
#'
#' @param lis list of indices
#' @param idx vector of indices
#' @param adj adjacency matrix
#'
#' @return one element of \code{lis}
#' @export
select_clique <- function(lis, idx, adj){
  stopifnot(max(idx) <= ncol(adj), ncol(adj) == nrow(adj))
  stopifnot(all(sapply(lis, max) <= ncol(adj)))
  
  if(length(lis) > 1){
    vec <- sapply(lis, function(x){
      length(which(idx %in% x))
    })
    
    if(length(which(vec == max(vec))) > 1){
      idx <- which.max(sapply(lis, function(x){sum(adj[x])}))
      lis[[idx]]
    } else {
      lis[[which(vec == max(vec))]]
    }
  } else {
    lis[[1]]
  }
}


#' Determine if an adjacency matrix passes a threshold
#'
#' Compares the number of non-zeros in an adjacency matrix to a threshold
#'
#' @param adj 0-1 symmetrical matrix
#' @param threshold Binomial probability
#'
#' @return logical
.pass_threshold <- function(adj, threshold){
  stopifnot(is.matrix(adj))
  
  if(length(adj) == 1) return(TRUE)
  
  n <- ncol(adj); n2 <- n*(n-1)/2
  num_edges <- sum(adj[upper.tri(adj)])
  cutoff <- n2*threshold
  if(num_edges > cutoff) return(TRUE) else return(FALSE)
}

.convert_string_to_pair <- function(string){
  vec <- as.numeric(strsplit(string, split = "-")[[1]])
  stopifnot(vec[1] < vec[2])
  vec
}

.convert_pair_to_string <- function(vec, pair = T){
  if(pair) stopifnot(length(vec) == 2)
  stopifnot(is.numeric(vec), !is.matrix(vec), length(unique(vec)) == length(vec))
  vec <- sort(vec)
  paste0(vec, collapse = "-")
}

.convert_indices_to_binary <- function(idx, max){
  stopifnot(all(idx <= max), all(idx > 0), all(idx %% 1 == 0))
  vec <- rep(0, max)
  vec[idx] <- 1
  paste0(vec, collapse = "")
}

.convert_binary_to_indices <- function(string){
  vec <- as.numeric(strsplit(string, split = "")[[1]])
  which(vec == 1)
}


#' Initialize queue
#'
#' Make a \code{queue} object and put all (len choose 2) strings into it
#'
#' @param len integer
#'
#' @return a \code{queue} object
.initialize_queue <- function(len){
  stopifnot(len > 1)
  
  queue <- dequer::queue()
  combn_mat <- utils::combn(len, 2)
  for(i in 1:ncol(combn_mat)){
    string <- .convert_pair_to_string(combn_mat[,i])
    dequer::pushback(queue, string)
  }
  
  queue
}

#' Initialize hash table of nodes
#'
#' Given a list of indices, make a node for each element and put it into a
#' hash table
#'
#' @param lis list of indices
#'
#' @return a \code{hash} object of \code{node} objects
.initialize_children <- function(lis){
  has <- hash::hash()
  for(i in 1:length(lis)){
    node <- .node(lis[[i]], NA)
    has[[as.character(i)]] <- node
  }
  
  has
}

#' Initialize hash table of indices
#'
#' @param lis list of indices
#' @param d maximum value in \code{elements}
#'
#' @return a \code{hash} object of booleans
.initialize_unique <- function(lis, d){
  has <- hash::hash()
  for(i in 1:length(lis)){
    string <- .convert_indices_to_binary(lis[[i]], max = d)
    has[[string]] <- TRUE
  }
  
  has
}

#' Create a \code{node} object
#'
#' A \code{node} object contains \code{elements} (the indices) and
#' \code{children} (the assigned number of the two node names that this
#' current node is made from). Note that this node name is not stored in this
#' object, but in how it is referenced (in a hash table, for instance).
#'
#' \code{children} can be either 2 positive integers (in which case the first number
#' must be small than the second) or \code{NA}.
#'
#' @param elements vector
#' @param children vector
#'
#' @return a \code{node} object
.node <- function(elements, children){
  stopifnot(is.numeric(elements), all(elements %% 1 == 0), all(elements > 0))
  stopifnot(is.na(children) || (length(children) == 2 & children[1] < children[2] &
                                  all(children %% 1 == 0) & all(children > 0)))
  
  structure(list(elements = elements, children = children), class = "node")
}

#' Check superset
#'
#' Return \code{TRUE} if either \code{elements1} or \code{elements2} is a
#' subset of another
#'
#' @param elements1 vector of indices
#' @param elements2 vector of indices
#'
#' @return boolean
.check_superset <- function(elements1, elements2){
  all(elements1 %in% elements2) | all(elements2 %in% elements1)
}

#' Check if an attempted set formed by two sets of indices has already been tried
#'
#' Function returns \code{FALSE} if the set is unique (meaning, it has not been
#' tried before), return \code{TRUE} if the set is not unique (meaning it has been
#' tried before).
#'
#' @param has hash table of booleans with strings of indicies as keys
#' @param elements1 vector of indices
#' @param elements2 vector of indices
#' @param d maximum value in \code{elements}
#'
#' @return boolean
.check_nonunique <- function(has, elements1, elements2, d){
  vec <- sort(unique(c(elements1, elements2)))
  
  string <- .convert_indices_to_binary(vec, max = d)
  val <- has[[string]]
  !is.null(val)
}

#' Check pairs of children
#'
#' Given two vectors, \code{child1} and \code{child2}, and a hash table of
#' booleans (\code{has}), check to see if any/all of the elements in \code{child1} with
#' any of the elements in \code{child2} are \code{TRUE} or \code{FALSE}.
#'
#' \code{child1} and \code{child2} may be \code{NA}, in which case this function
#' always returns a \code{TRUE}. Otherwise, \code{child1} and \code{child2}
#' can be a vector of positive integers with length 1 or 2.
#'
#' When mode "all" is used, a \code{TRUE} is returned only if all pairs between
#' \code{child1} and \code{child2} are \code{TRUE}. Mode "or" is when only one
#' pair needs to be \code{TRUE} to output \code{TRUE}.
#'
#' \code{null_alarm} is a boolean. If \code{TRUE}, this function will throw
#' an error if a pair is attempted to be checked, but it is not recorded in
#' \code{has}.
#'
#' @param child1 vector
#' @param child2 vector
#' @param has a \code{hash} object
#' @param mode string
#' @param null_alarm boolean
#'
#' @return boolean
.check_pairs <- function(child1, child2, has, mode = "all", null_alarm = F){
  stopifnot(is.numeric(child1) | is.na(child1))
  stopifnot(is.numeric(child2) | is.na(child2))
  stopifnot(length(child1) <= 2, length(child1) >= 1)
  stopifnot(length(child2) <= 2, length(child2) >= 1)
  
  mat <- expand.grid(child1, child2)
  
  bool_vec <- apply(mat, 1, function(x){
    if(any(is.na(x))) return(TRUE)
    if(length(unique(x)) == 1) return(TRUE)
    string <- .convert_pair_to_string(x)
    res <- has[[string]]
    if(any(is.null(res))) {
      if(null_alarm) stop("NULL found in .check_pairs") else return(TRUE)
    }
    
    stopifnot(is.logical(res), length(res) == 1)
    res
  })
  
  if(mode == "all") {
    all(bool_vec)
  } else {
    any(bool_vec)
  }
}

#' Add new elements to a queue
#'
#' This function pushes the strings formed by any integer less than \code{new_len}
#' and \code{new_len} onto the queue.
#'
#' Recall, since \code{queue} is using \code{dequer}'s implementation, it is
#' passed by reference.
#'
#' @param queue \code{queue} object
#' @param new_len positive integer
#' @param elements indices of the new node to be added
#' @param children vector
#' @param hash_children \code{hash} object
#' @param hash_unique \code{hash} object
#' @param hash_history \code{hash} object
#' @param d maximum value in \code{elements}
#' @param mode string
#'
#' @return a vector, stating how many things were attempted to be pushed and how many things were actually pushed
.add_to_queue <- function(queue, new_len, elements, children, hash_children, hash_unique,
                          hash_history, d, mode = "all"){
  stopifnot(length(new_len) == 1, is.numeric(new_len), new_len > 1)
  counter <- 0
  
  for(i in 1:(new_len - 1)){
    string <- .convert_pair_to_string(c(i, new_len))
    
    elements_other <- hash_children[[as.character(i)]]$elements
    bool <- .check_superset(elements, elements_other)
    if(bool) {
      hash_history[[string]] <- TRUE
      next()
    } else {
      bool <- .check_nonunique(hash_unique, elements, elements_other, d)
      if(bool) {
        hash_history[[string]] <- TRUE
        next()
      }
    }
    
    #list children of two pairs, check all pairs
    children2 <- hash_children[[as.character(i)]]$children
    bool <- .check_pairs(children, children2, hash_history, mode = mode, null_alarm = F)
    if(!bool) {
      hash_history[[string]] <- FALSE; next()
    }
    
    dequer::pushback(queue, string)
    counter <- counter + 1
  }
  
  c(new_len, counter)
}

#' Checks if all nodes in a hash table are unique
#'
#' @param hash a \code{hash} object
#' @param len maximum index
#'
#' @return boolean
.check_node_uniqeness <- function(hash, len){
  val <- hash::values(hash)
  stopifnot(all(rownames(val) == c("elements", "children")))
  
  mat <- sapply(1:ncol(val), function(x){
    vec <- rep(0, len)
    vec[val[1,x][[1]]] <- 1
    vec
  })
  
  mat <- t(mat)
  dis_mat <- stats::dist(mat)
  
  !any(dis_mat == 0)
}

#' Prune the list of cliques
#'
#' @param adj 0-1 symmetrical matrix
#' @param clique_list list of vectors
#' @param target_idx vector
#' @param threshold target density of quasi-clique
#'
#' @return list of vectors
.prune_clique <- function(adj, clique_list, target_idx, threshold = 0.9){
  vec <- sapply(clique_list, function(x){
    idx <- unique(c(x, target_idx))
    .pass_threshold(adj[idx,idx], threshold)
  })
  
  if(!any(vec)) stop("Pruning maximal cliques resulted in no sets")
  clique_list <- clique_list[which(vec)]
  
  lapply(clique_list, function(x){
    sort(unique(c(x, target_idx)))
  })
}

#from rje package, original function "powerSet"
.powerset <- function(x){
  out <- list(x[c()])
  if (length(x) == 1)
    return(c(out, list(x)))
  for (i in seq_along(x)) {
    out <- c(out, lapply(out, function(y) c(y, x[i])))
  }
  out
}

#' Remove sparsely connected nodes
#'
#' @param clique_list list of vectors
#' @param adj adjacency graph
#' @param target_idx vector of indices from 1 to \code{ncol(adj)} to force vectors in \code{clique_list} to retain
#' @param threshold number between 0 and 1
#'
#' @return list of vectors
.prune_nodes <- function(clique_list, adj, target_idx, threshold = 0.5){
  for(i in 1:length(clique_list)){
    if(length(clique_list[[i]]) == 0) next()
    adj_small <- adj[clique_list[[i]], clique_list[[i]], drop = F]
    vec <- rowSums(adj_small)
    
    if(!all(is.na(target_idx))){
      passed_idx <- sort(unique(c(which(vec >= threshold*nrow(adj_small)), which(clique_list[[i]] %in% target_idx))))
    } else {
      passed_idx <- which(vec >= threshold*nrow(adj_small))
    }
    
    clique_list[[i]] <- clique_list[[i]][passed_idx]
  }
  
  clique_list
}