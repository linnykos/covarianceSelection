context("Test extractor")

## .split_name is correct

test_that(".split_name works", {
  vec <- "HSB121.MD.R.CEL"
  res <- .split_name(vec)
  expect_true(all(res == c("HSB121", "MDCBC", "8")))
})

test_that(".split_name can handle the region itself", {
  vec <- "HSB100.PFC.6"
  res <- .split_name(vec)
  expect_true(all(res == c("HSB100", "PFC", "6")))
})

## extractor is correct

test_that("extractor works", {
  dat <- as.data.frame(matrix(rnorm(62), 31, 2))
  rownames(dat) <- c("HSB136.A1C.R.CEL",
                     "HSB136.AMY.L.CEL", "HSB136.AMY.R.CEL",
                     "HSB136.CBC.L.CEL", "HSB136.CBC.R.CEL",
                     "HSB136.DFC.L.CEL", "HSB136.DFC.R.CEL",
                     "HSB136.HIP.L.CEL", "HSB136.HIP.R.CEL",
                     "HSB136.IPC.L.CEL", "HSB136.IPC.R.CEL",
                     "HSB136.ITC.L.CEL", "HSB136.ITC.R.CEL",
                     "HSB136.M1C.L.CEL", "HSB136.M1C.R.CEL",
                     "HSB136.MD.L.CEL", "HSB136.MD.R.CEL",
                     "HSB136.MFC.L.CEL", "HSB136.MFC.R.CEL",
                     "HSB136.OFC.L.CEL", "HSB136.OFC.R.CEL",
                     "HSB136.S1C.L.CEL", "HSB136.S1C.R.CEL",
                     "HSB136.STC.L.CEL", "HSB136.STC.R.CEL",
                     "HSB136.STR.L.CEL", "HSB136.STR.R.CEL",
                     "HSB136.V1C.L.CEL", "HSB136.V1C.R.CEL",
                     "HSB136.VFC.L.CEL", "HSB136.VFC.R.CEL")
  res <- extractor(dat)

  expect_true(is.list(res))
  expect_true(length(res) == 4)
  expect_true(all(sapply(res, ncol) == 2))
  expect_true(all(sort(sapply(res, nrow)) == sort(c(12, 9, 6, 4))))
})
