context("Test isCheck, type checkers")

## test .is.listOfMatrix

test_that("it errors on non-lists", {
  expect_error(.is.listOfMatrix(1:5, "lis"))
})

test_that("it errors when the list has non-matrices", {
  lis <- vector("list", 2)
  lis[[1]] <- matrix(1:25, 5, 5)
  lis[[2]] <- 1:5
  expect_error(.is.listOfMatrix(lis, "lis"))
})

test_that("it works on list of matrices", {
  lis <- vector("list", 2)
  lis[[1]] <- matrix(1:25, 5, 5)
  lis[[2]] <- matrix(1:12, 4, 3)
  expect_true(.is.listOfMatrix(lis, "lis"))
})

##############################################

## test .is.listofNumeric

test_that("it errors on non-lists", {
  expect_error(.is.listofNumeric(1:5, "lis"))
})

test_that("it errors when the list has non-numerics", {
  lis <- vector("list", 2)
  lis[[1]] <- matrix(1:25, 5, 5)
  lis[[2]] <- 1:5
  expect_error(.is.listofNumeric(lis, "lis"))
})

test_that("it works on list of matrices", {
  lis <- vector("list", 2)
  lis[[1]] <- 1:5
  lis[[2]] <- 1:10
  expect_true(.is.listofNumeric(lis, "lis"))
})

test_that("it errors when there are strings", {
  lis <- vector("list", 2)
  lis[[1]] <- c("a","b","c")
  lis[[2]] <- 1:5
  expect_error(.is.listofNumeric(lis, "lis"))
})

############################################3

## test .is.nonNegInteger

test_that("it errors on NAs with numerics", {
  expect_error(.is.nonNegInteger(c(3, NA, 5), "vec"))
  expect_error(.is.nonNegInteger(c(NA, 3, 5), "vec"))
})

test_that("it works on numerics",{
  expect_true(.is.nonNegInteger(1:10, "vec"))
})

test_that("it errors on negatives", {
  expect_error(.is.nonNegInteger(-1:5, "vec"))
})

test_that("it errors on decimals", {
  expect_error(.is.nonNegInteger(1:10+.5, "vec"))
  expect_error(.is.nonNegInteger(c(1,6,4,2.4,2), "vec"))
})

test_that("it errors on duplicates", {
  expect_error(.is.nonNegInteger(c(1,5,5), "vec"))
})
