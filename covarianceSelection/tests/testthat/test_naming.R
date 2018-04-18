context("Test naming")
#
# ## convert_genenames is correct
#
# test_that("convert_genenames works", {
#   library(biomaRt)
#   vec <- c("KDM6B", "NCKAP1", "RANBP17", "WDFY3", "DSCAM",
#            "FOXP1", "KDM5B", "MED13L", "PHF2", "SPAST", "CHD8", "DYRK1A",
#            "ANK2" , "ARID1B", "ADNP", "POGZ")
#   res <- convert_genenames(vec)
#
#   expect_true(length(res) == length(vec))
#   expect_true(all(sapply(res, function(x){substr(x, 1, 4)}) == "ENSG"))
# })
#
# test_that("convert_genenames is reflective", {
#   library(biomaRt)
#   vec <- c("KDM6B", "NCKAP1", "RANBP17", "WDFY3", "DSCAM",
#            "FOXP1", "KDM5B", "MED13L", "PHF2", "SPAST", "CHD8", "DYRK1A",
#            "ANK2" , "ARID1B", "ADNP", "POGZ")
#   res <- convert_genenames(vec)
#   res2 <- convert_genenames(res, from = "ensembl", to = "symbol")
#
#   expect_true(all(vec == res2))
# })
#
# test_that("convert_genenames can output NA", {
#   library(biomaRt)
#   vec <- c("asdf", "KDM6B", "NCKAP1", "RANBP17", "WDFY3", "DSCAM",
#            "FOXP1", "KDM5B", "MED13L", "PHF2", "SPAST", "CHD8", "DYRK1A",
#            "ANK2" , "ARID1B", "ADNP", "asdfasdf", "POGZ")
#   res <- convert_genenames(vec)
#
#   expect_true(all(is.na(res[c(1,17)])))
#   expect_true(all(!is.na(res[-c(1,17)])))
# })

###########################

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
