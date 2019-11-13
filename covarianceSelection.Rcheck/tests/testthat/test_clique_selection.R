context("Test clique selection")

## .convert_string_to_pair is correct

test_that(".convert_string_to_pair works", {
  res <- .convert_string_to_pair("4-10")
  expect_true(all(res == c(4,10)))
})

test_that(".convert_string_to_pair is reversible", {
  string <- "67-123"
  string2 <- .convert_pair_to_string(.convert_string_to_pair(string))
  
  expect_true(string == string2)
})


##########################

## .convert_pair_to_string is correct

test_that(".convert_pair_to_string works", {
  res <- .convert_pair_to_string(c(4,10))
  expect_true(res == "4-10")
})

test_that(".convert_pair_to_string is reversible", {
  vec <- c(27,67)
  vec2 <- .convert_string_to_pair(.convert_pair_to_string(vec))
  
  expect_true(all(vec == vec2))
})

##########################

## .initialize_queue is correct

test_that(".initialize_queue works", {
  res <- .initialize_queue(10)
  
  expect_true(dequer:::length.queue(res) == 10*9/2)
  
  vec <- rep(NA, 45)
  for(i in 1:45){
    vec[i] <- dequer::pop(res)
  }
  
  expect_true(dequer:::length.queue(res) == 0)
  expect_true(length(unique(vec)) == length(vec))
})

########################

## .node is correct

test_that(".node works", {
  res <- .node(1:5, c(4,8))
  expect_true(class(res) == "node")
  expect_true(all(res$elements == 1:5))
  expect_true(all(res$children == c(4,8)))
})

test_that(".node can take no children", {
  res <- .node(1:5, NA)
  
  expect_true(class(res) == "node")
  expect_true(is.na(res$children))
})

 ########################

## .initialize_children is correct

test_that(".initialize_children works", {
  load("../assets/clique_selection1.RData")
  
  len <- 25
  edges <- combn(25, 2)[,lis[[1]]]
  g <- igraph::graph.empty(n = len, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  
  lis <- lapply(igraph::maximal.cliques(g), as.numeric)
  
  res <- .initialize_children(lis)
  
  expect_true(length(res) == length(lis))
  class_vec <- sapply(1:length(lis), function(x){
    class(res[[as.character(x)]])
  })
  expect_true(all(class_vec == "node"))
})

############################

## .check_superset is correct

test_that(".check_superset works", {
  load("../assets/clique_selection1.RData")
  
  res <- .check_superset(c(1:10), seq(1,9,by=2))
  
  expect_true(length(res) == 1)
  expect_true(is.logical(res))
})

##########################


## .check_pairs is correct

test_that(".check_pairs works", {
  hash_history <- hash::hash()
  res <- .check_pairs(NA, NA, hash_history)
  
  expect_true(length(res) == 1)
  expect_true(is.logical(res))
})

test_that(".check_pairs works for pairs of children, all tested", {
  hash_history <- hash::hash()
  hash_history[["1-3"]] <- TRUE
  hash_history[["1-4"]] <- TRUE
  hash_history[["2-3"]] <- TRUE
  hash_history[["2-4"]] <- TRUE
  
  res <- .check_pairs(c(1,2), c(3,4), hash_history)
  
  expect_true(res)
})

test_that(".check_pairs works for pairs of children, some tested", {
  hash_history <- hash::hash()
  hash_history[["1-3"]] <- TRUE
  hash_history[["1-4"]] <- TRUE
  hash_history[["2-3"]] <- TRUE
  
  res <- .check_pairs(c(1,2), c(3,4), hash_history)
  
  expect_true(res)
})

test_that(".check_pairs works for known failed pairs, using mode ALL", {
  hash_history <- hash::hash()
  hash_history[["1-3"]] <- TRUE
  hash_history[["1-4"]] <- TRUE
  hash_history[["2-3"]] <- FALSE
  
  res <- .check_pairs(c(1,2), c(3,4), hash_history)
  
  expect_true(!res)
})

test_that(".check_pairs works for known failed pairs, using mode OR", {
  hash_history <- hash::hash()
  hash_history[["1-3"]] <- TRUE
  hash_history[["1-4"]] <- TRUE
  hash_history[["2-3"]] <- FALSE
  
  res <- .check_pairs(c(1,2), c(3,4), hash_history, mode = "or")
  
  expect_true(res)
})

test_that(".check_pairs can flag NA's", {
  hash_history <- hash::hash()
  hash_history[["1-3"]] <- TRUE
  hash_history[["1-4"]] <- TRUE
  hash_history[["2-3"]] <- TRUE
  
  expect_error(.check_pairs(c(1,2), c(3,4), hash_history, null_alarm = T))
})

############

## .add_to_queue is correct

test_that(".add_to_queue works", {
  n <- 10
  combn_mat <- combn(n,2)
  
  set.seed(10)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  
  queue <- .initialize_queue(n)
  hash_history <- hash::hash()
  clique_list <- lapply(igraph::maximal.cliques(g), function(x){sort(as.numeric(x))})
  hash_children <- .initialize_children(clique_list)
  hash_unique <- .initialize_unique(clique_list, n)
  
  len <- length(clique_list)
  .add_to_queue(queue, len+1, sort(unique(unlist(clique_list[1:2]))), c(1,2),
                hash_children, hash_unique, hash_history, n)
  
  expect_true(length(queue) > len)
})

################

## clique_selection is correct

test_that("clique_selection works", {
  set.seed(10)
  combn_mat <- combn(10,2)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
  
  g <- igraph::graph.empty(n = 10, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  
  res <- clique_selection(g)
  
  expect_true(is.numeric(res[[1]]))
  expect_true(!is.matrix(res[[1]]))
  expect_true(all(res[[1]] %% 1 == 0))
  expect_true(all(res[[1]] >= 1))
  expect_true(is.list(res))
})

test_that("clique_selection works with the heuristic", {
  set.seed(10)
  combn_mat <- combn(10,2)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
  
  g <- igraph::graph.empty(n = 10, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  
  res <- clique_selection(g, target_idx = c(1:5))
  
  expect_true(is.list(res))
})

test_that("clique_selection does not crash", {
  trials <- 20
  combn_mat <- combn(10,2)
  
  bool_vec <- rep(NA, trials)
  for(i in 1:trials){
    set.seed(i)
    edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
    
    g <- igraph::graph.empty(n = 10, directed = F)
    g <- igraph::add_edges(g, edges = edges)
    
    res <- clique_selection(g)[[1]]
    
    bool_vec[i] <- is.numeric(res)
  }
  
  expect_true(all(bool_vec))
})

test_that("clique_selection works on full cliques", {
  g <- igraph::graph.empty(n = 10, directed = F)
  g <- igraph::add_edges(g, edges = combn(10,2))
  res <- clique_selection(g)[[1]]
  
  expect_true(all(res == 1:10))
})

test_that("clique_selection works with on empty graphs", {
  n <- 20
  g <- igraph::graph.empty(n = n, directed = F)
  res <- clique_selection(g)
  
  expect_true(length(res[[1]]) == 0)
})

test_that("clique_selection can find large clique", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), 
                    tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- clique_selection(g)
  
  expect_true(all(res[[1]] == (n/2+1):n))
})

test_that("clique_selection all pass threshold", {
  trials <- 20
  combn_mat <- combn(10,2)
  
  bool_vec <- rep(NA, trials)
  for(i in 1:trials){
    set.seed(20*i)
    edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
    
    g <- igraph::graph.empty(n = 10, directed = F)
    g <- igraph::add_edges(g, edges = edges)
    
    res <- clique_selection(g, threshold = 0.95)[[1]]
    
    adj <- as.matrix(igraph::as_adjacency_matrix(g))
    bool_vec[i] <- .pass_threshold(adj[res, res, drop = F], threshold = 0.95)
  }
  
  expect_true(all(bool_vec))
})


test_that("clique_selection does not repeat nodes", {
  trials <- 20
  combn_mat <- combn(10,2)
  
  bool_vec <- rep(NA, trials)
  for(i in 1:trials){
    set.seed(20*i)
    edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
    
    g <- igraph::graph.empty(n = 10, directed = F)
    g <- igraph::add_edges(g, edges = edges)
    
    res <- clique_selection(g, threshold = 0.95)[[1]]
    
    adj <- as.matrix(igraph::as_adjacency_matrix(g))
    bool_vec[i] <- .pass_threshold(adj[res, res, drop = F], threshold = 0.95)
  }
  
  expect_true(all(bool_vec))
})

test_that("clique_selection does not have uniqueness problems", {
  n <- 13
  combn_mat <- combn(n,2)
  
  set.seed(10)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  res <- clique_selection(g)
  
  expect_true(is.list(res))
})

test_that("clique_selection does not suffer overflow problem", {
  load("../assets/clique_selection2.RData")
  
  g <- igraph::graph.empty(n = 25, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  res <- clique_selection(g)
  
  expect_true(is.list(res))
})

test_that("clique_selection gives the proper output", {
  load("../assets/clique_selection3.RData")
  
  num_partition <- 25
  edges <- combn(num_partition, 2)
  g <- igraph::graph.empty(n = num_partition, directed = F)
  g <- igraph::add_edges(g, edges = edges[, indices_list[[8]]])
  res <- clique_selection(g, threshold = threshold,
                          mode = "or", verbose = F,
                          time_limit = 180)
  res <- select_clique(res, 1:15, g)
  
  idx <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
           15, 16, 18, 19, 20, 21, 22, 23, 25)
  expect_true(all(sort(res) == idx))
  adj <- as.matrix(igraph::as_adj(g))
  expect_true(.pass_threshold(adj[idx, idx], 0.95))
})

test_that("clique_selection is not random", {
  load("../assets/clique_selection3.RData")
  
  trials <- 100
  mat <- sapply(1:trials, function(x){
    num_partition <- 25
    edges <- combn(num_partition, 2)
    g <- igraph::graph.empty(n = num_partition, directed = F)
    g <- igraph::add_edges(g, edges = edges[, indices_list[[8]]])
    res <- clique_selection(g, threshold = threshold,
                            mode = "or", verbose = F,
                            time_limit = 180)
    select_clique(res, 1:15, g)
  })
  
  bool_vec <- apply(mat, 1, function(x){length(unique(x)) == 1})
  expect_true(all(bool_vec))
})


#######################

## .pass_threshold is correct

test_that(".pass_threshold works", {
  set.seed(10)
  adj_mat <- matrix(sample(c(0,1),100, replace = T), 10, 10)
  adj_mat <- pmax(adj_mat, t(adj_mat))
  
  res <- .pass_threshold(adj_mat, 0.5)
  
  expect_true(is.logical(res))
  expect_true(length(res) == 1)
})

test_that(".pass_threshold can accept", {
  adj_mat <- matrix(1, 10, 10)
  res <- .pass_threshold(adj_mat, 0.95)
  expect_true(res)
})

test_that(".pass_threshold can reject", {
  set.seed(10)
  adj_mat <- matrix(sample(c(0,1),100, replace = T), 10, 10)
  adj_mat <- adj_mat * t(adj_mat)
  
  res <- .pass_threshold(adj_mat, 0.95)
  
  expect_true(!res)
})


###############################

## select_clique is correct

test_that("select_clique works", {
  set.seed(10)
  combn_mat <- combn(11,2)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.9*ncol(combn_mat)))]
  
  g <- igraph::graph.empty(n = 11, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  
  clique_list <- clique_selection(g)
  idx <- select_clique(clique_list, 1:5, g)
  
  expect_true(is.numeric(idx))
  expect_true(!is.matrix(idx))
  expect_true(!is.list(idx))
})

#################################

## .prune_clique is correct

test_that(".prune_clique works", {
  set.seed(10)
  n <- 20
  combn_mat <- combn(n,2)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.7*ncol(combn_mat)))]
  
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  adj <- as.matrix(igraph::as_adjacency_matrix(g))
  
  clique_list <- lapply(igraph::maximal.cliques(g), function(x){sort(as.numeric(x))})
  
  res <- .prune_clique(adj, clique_list, target_idx = 1:10, threshold = 0.5)
  
  expect_true(is.list(res))
})

test_that(".prune_clique returns a strict subset, but with target idx appended", {
  set.seed(10)
  n <- 20
  combn_mat <- combn(n,2)
  edges <- combn_mat[,sample(1:ncol(combn_mat), floor(0.7*ncol(combn_mat)))]
  
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  adj <- as.matrix(igraph::as_adjacency_matrix(g))
  
  clique_list <- lapply(igraph::maximal.cliques(g), function(x){sort(as.numeric(x))})
  
  res <- .prune_clique(adj, clique_list, target_idx = 1:10, threshold = 0.5)
  
  clique_string <- sapply(clique_list, function(x){
    paste0(sort(unique(c(x, 1:10))), collapse = "-")
  })
  res_string <- sapply(res, paste0, collapse = "-")
  bool <- res_string %in% clique_string
  
  expect_true(all(bool))
})
