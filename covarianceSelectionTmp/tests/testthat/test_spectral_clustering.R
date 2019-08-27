context("Test spectral clustering")

## spectral_selection works

test_that("spectral_selection works", {
  load("../assets/clique_selection3.RData")

  num_partition <- 25
  edges <- combn(num_partition, 2)
  g <- igraph::graph.empty(n = num_partition, directed = F)
  g <- igraph::add_edges(g, edges = edges[, indices_list[[8]]])

  set.seed(10)
  res <- spectral_selection(g, K = 3, target_idx <- c(1:3, 16, 21))

  expect_true(length(res) <= num_partition)
})
