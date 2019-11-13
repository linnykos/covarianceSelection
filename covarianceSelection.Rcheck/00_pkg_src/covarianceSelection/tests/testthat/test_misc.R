context("Test miscellaenous functions")

## matching is correct

test_that("matching works", {
  set.seed(10)
  name1 <- letters[sample(1:26)]
  name2 <- letters[sample(1:26)]
  
  res <- matching(name1, name2)
  
  expect_true(length(res) == 26)
  expect_true(all(res > 0))
  expect_true(all(res <= 26))
  expect_true(all(res %% 1 == 0))
  expect_true(length(res) == length(unique(res)))
})

test_that("matching is correct, simple example", {
  set.seed(10)
  vec1 <- rnorm(10)
  vec2 <- vec1[c(1,3,5,2,6,7,8,4,9,10)]
  
  res <- matching(vec1, vec2)
  
  expect_true(all(vec1 == vec2[res]))
})
