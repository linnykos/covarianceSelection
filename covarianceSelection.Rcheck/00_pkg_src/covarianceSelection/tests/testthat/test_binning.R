context("Test binning")

## binning is correct

test_that("binning works", {
  vec <- c("HSB154.PFC.5", "HSB97.PFC.5", "HSB100.PFC.6", "HSB100.VIIAS.6")
  res <- binning(vec)

  expect_true(all(dim(res) == c(4,4)))
  expect_true(all(res >= 0))
  expect_true(all(res <= length(vec)))
  expect_true(all(as.numeric(res)%%1 == 0))
})

test_that("binning gets the right results", {
  vec <- c("HSB154.PFC.5", "HSB97.PFC.5", "HSB100.PFC.6", "HSB100.VIIAS.6")
  res <- binning(vec)

  answer <- c(0,0,0,0, 2,0,0,0, 1,1,0,0, 0,0,0,0)
  expect_true(all(as.numeric(res) == answer))
})

