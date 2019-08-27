context("Test graph functions")

## compute_scale_free is correct

test_that("compute_scale_free works", {
  set.seed(10)
  d <- 50
  dat <- huge::huge.generator(n = 50, d = d, graph = "hub", g = 4, verbose = F)
  adj <- as.matrix(dat$theta)

  res <- compute_scale_free(adj)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

##################

## enumerate_edges is correct

test_that("enumerate_edges works", {
  adj <- matrix(c(1, 0, 0, 1,
                  0, 1, 1, 0,
                  0, 1, 1, 1,
                  1, 0, 1, 1), 4, 4)
  res <- enumerate_edges(adj)

  expect_true(is.matrix(res))
  expect_true(ncol(res) == 2)
  expect_true(nrow(res) == 3)
})
