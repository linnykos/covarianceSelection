context("Test goodness of fit")

## .partition_data is correct

test_that(".partition_data works", {
  set.seed(1)
  res <- .partition_data(10, 5)

  expect_true(is.list(res))
  expect_true(all(sapply(res, is.list)))
  expect_true(length(res) == 5)
  expect_true(all(sapply(res, length) == 2))
  expect_true(all(sapply(res, function(x){
    length(unique(unlist(x)))
  }) == 10))
  expect_true(all(sapply(res, function(x){
    length(unlist(x))
  }) == 10))
})

test_that(".partition_data can return all partitions", {
  set.seed(1)
  n <- 5
  res <- .partition_data(n, 20)

  expect_true(length(res) == 2^4-1)

  bool_vec <- sapply(1:20, function(x){
    set.seed(x)
    y <- ceiling(runif(1, 0, n-1))
    left <- sample(1:5)[y]

    vec <- sapply(res, function(y){
      all(y[[1]] == left) | all(y[[2]] == left)
    })

    sum(vec) == 1
  })

  expect_true(all(bool_vec))
})

test_that(".partition_data returns only unique partitions, even when exceeding max", {
  set.seed(1)
  n <- 5
  res <- .partition_data(n, 20)

  vec <- sapply(res, function(x){
    if(all(x[[1]] == n)) return(0)
    sum(sapply(x[[1]], function(y){2^(y-1)}))
  })

  expect_true(length(unique(vec)) == 2^(n-1)-1)
  expect_true(all(sort(vec) == c(1:(2^(n-1)-1))))
})

test_that(".partition_data returns only unique partitions", {
  set.seed(1)
  n <- 10
  res <- .partition_data(n, 100)

  vec <- sapply(res, function(x){
    if(all(x[[1]] == n)) return(0)
    sum(sapply(x[[1]], function(y){2^(y-1)}))
  })

  expect_true(length(unique(vec)) == 100)
})

test_that(".partition_data does not crash for large n's", {
  set.seed(1)
  n <- 50
  res <- .partition_data(n, 100)

  expect_true(is.list(res))
})


################

## goodness_of_fit is correct

test_that("goodness_of_fit works", {
  set.seed(1)
  dat_list <- lapply(1:10, function(x){
    MASS::mvrnorm(50, mu = rep(0,4), Sigma = diag(4))
  })

  res <- goodness_of_fit(dat_list, 10, trials = 15)

  expect_true(length(res) == 10)
  expect_true(all(res >= 0))
  expect_true(all(res <= 1))
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})
