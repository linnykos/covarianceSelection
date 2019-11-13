context("Test stepdown accelerated")

## .construct_graph is correct

test_that(".construct_graph works", {
  g <- .construct_graph(5)

  expect_true(class(g) == "igraph")
})

test_that(".construct_graph can construct a subset graph", {
  edge_mat <- combn(10, 2)
  idx <- which(apply(edge_mat, 2, function(x){any(x == 2)}))
  edge_mat <- edge_mat[,-idx]

  g <- .construct_graph(edge_mat)

  expect_true(igraph::vcount(g) == 9)
  expect_true(igraph::ecount(g) == 9*8/2)
  expect_true(all(as.character(sort(as.numeric(igraph::V(g)$name))) == as.character(c(1,c(3:10)))))
})

test_that(".construct_graph sets the weight correctly", {
  edge_mat <- combn(10, 2)
  idx <- which(apply(edge_mat, 2, function(x){any(x == 2)}))
  edge_mat <- edge_mat[,-idx]

  g <- .construct_graph(edge_mat, value = Inf)

  expect_true(all(igraph::E(g)$weight == Inf))
})


test_that(".construct_graph works for incomplete graphs", {
  set.seed(10)

  edge_mat <- combn(10, 2)
  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat2 <- edge_mat[,idx]

  g <- .construct_graph(edge_mat2)

  bool_vec <- sapply(1:ncol(edge_mat), function(x){
    if(x %in% idx){
      vec <- range(edge_mat[,x])
      g[vec[1], vec[2]] == 1 & g[vec[2], vec[1]] == 1
    } else {
      vec <- range(edge_mat[,x])
      g[vec[1], vec[2]] == 0 & g[vec[2], vec[1]] == 0
    }
  })

  expect_true(all(bool_vec))
})


######

## .return_vertices is correct

test_that(".return_vertices works", {
  edge_mat <- combn(10, 2)
  idx <- which(apply(edge_mat, 2, function(x){any(x == 2)}))
  edge_mat <- edge_mat[,-idx]

  g <- .construct_graph(edge_mat)

  res <- .return_vertices(g)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(2, igraph::ecount(g))))
  expect_true(all(apply(res, 1, class) == "character"))
})

test_that(".return_vertices works for incomplete graphs", {
  set.seed(10)

  edge_mat <- combn(10, 2)
  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat2 <- edge_mat[,idx]

  g <- .construct_graph(edge_mat2)

  res <- .return_vertices(g)

  bool_vec <- sapply(1:ncol(edge_mat), function(x){
    if(x %in% idx){
      vec <- range(edge_mat[,x])
      length(intersect(which(res[1,] == vec[1]), which(res[2,] == vec[2]))) == 1
    } else {
      vec <- range(edge_mat[,x])
      length(intersect(which(res[1,] == vec[1]), which(res[2,] == vec[2]))) == 0
    }
  })

  expect_true(all(bool_vec))
})

#########

## .compute_difference_from_edges is correct

test_that(".compute_difference_from_edges works", {
  set.seed(10)
  n <- 10
  num_list <- lapply(1:n, function(x){
    matrix(rnorm(25), 5, 5)
  })

  edge_mat <- combn(10, 2)
  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat2 <- edge_mat[,idx]

  g <- .construct_graph(edge_mat2)
  vertex_mat <- .return_vertices(g)

  res <- .compute_difference_from_edges(num_list, vertex_mat)

  expect_true(length(res) == ncol(edge_mat2))
  expect_true(is.numeric(res))
})

#############

## .set_edge_weight is correct

test_that(".set_edge_weight works", {
  set.seed(10)
  n <- 10
  edge_mat <- combn(n, 2)
  g <- .construct_graph(edge_mat, value = Inf)

  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat2 <- edge_mat[,idx]
  g2 <- .construct_graph(edge_mat2)
  vertex_mat <- .return_vertices(g2)

  res <- .set_edge_weight(g, vertex_mat, 1:ncol(vertex_mat))

  expect_true(class(res) == "igraph")
})

test_that(".set_edge_weight sets the weights", {
  set.seed(20)
  n <- 10
  edge_mat <- combn(n, 2)
  g <- .construct_graph(edge_mat, value = Inf)

  expect_true(length(unique(igraph::E(g)$weight)) == 1)

  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat2 <- edge_mat[,idx]
  g2 <- .construct_graph(edge_mat2)
  vertex_mat <- .return_vertices(g2)

  res <- .set_edge_weight(g, vertex_mat, 1:ncol(vertex_mat))

  expect_true(length(unique(igraph::E(res)$weight)) > 1)
})

##########

## .compute_max_accelerated is correct

test_that(".compute_max_accelerated works", {
  set.seed(10)
  n <- 10
  num_list <- lapply(1:n, function(x){
    matrix(rnorm(25), 5, 5)
  })

  combn_mat <- combn(n, 2)

  res <- .compute_max_accelerated(num_list, combn_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".compute_max_accelerated can handle incomplete graphs", {
  set.seed(20)
  n <- 10
  num_list <- lapply(1:n, function(x){
    matrix(rnorm(25), 5, 5)
  })

  edge_mat <- combn(n, 2)
  idx <- sample(1:ncol(edge_mat), 20)
  edge_mat <- edge_mat[,idx]

  res <- .compute_max_accelerated(num_list, edge_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".compute_max_accelerated gets the same answer as .compute_all_test_stat_bootstrap", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    n <- 10
    num_list <- lapply(1:n, function(x){
      matrix(rnorm(25), 5, 5)
    })

    edge_mat <- combn(n, 2)
    idx <- sample(1:ncol(edge_mat), 20)
    edge_mat <- edge_mat[,idx]

    denom_list <- lapply(1:n, function(x){1})

    val1 <- max(abs(.compute_all_test_stat(num_list, denom_list, combn_mat = edge_mat, squared = F)))
    val2 <- .compute_max_accelerated(num_list, edge_mat)

    abs(val1 - val2) < 1e-6
  })

  expect_true(all(bool_vec))
})


