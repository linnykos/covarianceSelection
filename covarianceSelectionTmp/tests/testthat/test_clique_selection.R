context("Test clique selection")

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
  res <- clique_selection(g, num_pos = 1)[[1]]

  expect_true(all(res == 1:10))
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
  res <- clique_selection(g)
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
    res <- clique_selection(g, num_pos = 3, num_neg = 2)
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
