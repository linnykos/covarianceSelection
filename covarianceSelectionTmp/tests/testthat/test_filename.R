context("Test filename functions")

## filename_closure is correct

test_that("filename_closure works", {
  folder <- "asdf/"
  prefix <- "hi_"
  filename_func <- filename_closure(folder, prefix)

  res <- filename_func(10)

  expect_true(res == "asdf/hi_10.RData")
  expect_true(dir.exists(folder))

  unlink(folder, recursive = T)
  expect_true(!dir.exists(folder))
})
