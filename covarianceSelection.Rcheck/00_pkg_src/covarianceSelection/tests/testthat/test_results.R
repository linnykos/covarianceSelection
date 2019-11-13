context("Test results")

## report_results is correct

test_that("report_results works", {
  gene <- c("a", "b", "c", "d", "e")
  posterior <- stats::runif(5)
  pvalue <- stats::runif(5)
  Iupdate <- c(0,0,0,1,1)

  res <- report_results(gene, posterior, pvalue, Iupdate)
  expect_true(is.data.frame(res))
  expect_true(all(dim(res) == c(5,4)))
})

#############

## enrichment_test is correct

test_that("enrichment_test works", {
  vec1 <- as.character(1:5)
  vec2 <- as.character(seq(1,9,length.out = 5))
  all_vec <- as.character(1:15)

  res <- enrichment_test(vec1, vec2, all_vec)

  expect_true(is.list(res))
  expect_true(all(dim(res$contigency) == c(2,2)))
})
