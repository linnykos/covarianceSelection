context("Test pairwise test")

## .convert_combn2mat is correct

test_that(".convert_combn2mat works", {
  d <- 5
  mat <- matrix(1:25, d, d)
  comb <- combn(d, 2)
  vec <- apply(comb, 2, function(x){mat[x[2],x[1]]})
  vec2 <- mat[lower.tri(mat)]

  expect_true(all(vec == vec2))

  mat2 <- .convert_combn2mat(vec, 5)

  expect_true(all(mat[lower.tri(mat)] == mat2[lower.tri(mat2)]))
})
