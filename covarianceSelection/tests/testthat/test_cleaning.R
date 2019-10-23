context("Test cleaning")

## symbol_synonyms is correct

test_that("symbol_synonyms works", {
  vec <- c("ANKRD30BP2", "C21orf99", "CT85", "CTSP-1")
  res <- symbol_synonyms(vec, verbose = F)
  
  expect_true(length(unique(res)) == 1)
})

test_that("symbol_synonyms works with NAs", {
  vec <- c("ANKRD30BP2", "C21orf99", "CT85", "asdfasdfasdf", "CTSP-1")
  res <- symbol_synonyms(vec, verbose = F)
  
  expect_true(is.na(res[4]))
  expect_true(length(unique(res[c(1,2,3,5)])) == 1)
})

test_that("symbol_synonyms works with only NAs", {
  vec <- c("asdf", "asdfasdf")
  res <- symbol_synonyms(vec, verbose = F)
  
  expect_true(all(is.na(res)))
  expect_true(length(res) == 2)
})

###############


## test regroup

test_that("regroup properly partitions", {
  dat.list <- vector("list", 5)
  dat.list[[1]] <- matrix(1:25,5,5)
  dat.list[[2]] <- matrix(26:50,5,5)
  dat.list[[3]] <- matrix(51:75,5,5)
  dat.list[[4]] <- matrix(76:100,5,5)
  dat.list[[5]] <- matrix(101:125,5,5)
  
  idx.list <- vector("list", 2)
  idx.list[[1]] <- c(5,1,2)
  idx.list[[2]] <- c(3,4)
  
  res <- regroup(dat.list, idx.list)
  
  expect_true(length(res) == 2)
  expect_true(is.list(res))
  expect_true(all(sapply(res, is.matrix)))
  
  expect_true(all(dim(res[[1]]) == c(5,15)))
  expect_true(all(dim(res[[2]]) == c(5,10)))
  
  expect_true(all(res[[1]] == cbind(matrix(1:25,5,5), matrix(26:50,5,5),
                                    matrix(101:125,5,5))))
  expect_true(all(res[[2]] == cbind(matrix(51:75,5,5), matrix(76:100,5,5))))
})

## test .partition_renamed

test_that(".partition_renamed correctly converts", {
  idx.list <- vector("list", 4)
  idx.list[[1]] <- c(1,3)
  idx.list[[2]] <- 6
  idx.list[[3]] <- c(2,4,5)
  idx.list[[4]] <- 7:9
  
  vec <- .partition_renamed(idx.list)
  expect_true(length(vec) == 9)
  expect_true(all(vec == c(1,3,1,3,3,2,4,4,4)))
})

########################

## test .populate

test_that(".populate correctly populates", {
  idx.list <- vector("list", 4)
  idx.list[[1]] <- c(1,3)
  idx.list[[2]] <- 6
  idx.list[[3]] <- c(2,4,5)
  idx.list[[4]] <- 7:9
  
  ncol.vec <- 1:9
  
  vec <- .populate(ncol.vec, idx.list)
  expect_true(length(vec) == sum(ncol.vec))
  expect_true(all(vec == c(1, rep(3,2) , rep(1,3), rep(3,4), rep(3,5),
                           rep(2,6), rep(4, 7+8+9))))
})

############################

## test .regroup_check

test_that(".regroup_check errors if idx.list isn't a vector of ingeters", {
  expect_error(.regroup_check(10, c(1:5, "a")))
})

test_that(".regroup_check errors if idx.list is missing elements", {
  idx.list <- vector("list", 3)
  idx.list[[1]] <- c(1,3)
  idx.list[[2]] <- 6
  idx.list[[3]] <- c(2,4)
  
  expect_error(.regroup_check(6, idx.list))
})

test_that(".regroup_check errors if not consecutive", {
  idx.list <- vector("list", 3)
  idx.list[[1]] <- c(1,3)
  idx.list[[2]] <- 9
  idx.list[[3]] <- c(2,4,5)
  
  expect_error(.regroup_check(6, idx.list))
})

test_that(".regroup_check works as normal", {
  idx.list <- vector("list", 3)
  idx.list[[1]] <- c(1,3)
  idx.list[[2]] <- 6
  idx.list[[3]] <- c(2,4,5)
  
  expect_true(.regroup_check(6, idx.list))
})
