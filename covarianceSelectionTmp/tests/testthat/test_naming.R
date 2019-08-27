context("Test naming")

## symbol_synonyms is correct

test_that("symbol_synonyms works", {
  vec <- c("ANKRD30BP2", "C21orf99", "CT85", "CTSP-1")
  res <- symbol_synonyms(vec)

  expect_true(length(unique(res)) == 1)
})

test_that("symbol_synonyms works with NAs", {
  vec <- c("ANKRD30BP2", "C21orf99", "CT85", "asdfasdfasdf", "CTSP-1")
  res <- symbol_synonyms(vec)

  expect_true(is.na(res[4]))
  expect_true(length(unique(res[c(1,2,3,5)])) == 1)
})

test_that("symbol_synonyms works with only NAs", {
  vec <- c("asdf", "asdfasdf")
  res <- symbol_synonyms(vec)

  expect_true(all(is.na(res)))
  expect_true(length(res) == 2)
})
