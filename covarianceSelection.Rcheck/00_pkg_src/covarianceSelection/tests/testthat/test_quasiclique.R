context("Test quasiclique functions")

## tsourakakis_2014_approximate is correct
test_that("tsourakakis_2014_approximate works", {
  n <- 10
  g <- igraph::make_ring(n)
  res <- tsourakakis_2014_approximate(g)
  
  expect_true(is.numeric(res))
  expect_true(length(res) <= n)
  expect_true(length(unique(res)) == length(res))
  expect_true(all(res > 0))
  expect_true(all(res <= n))
  expect_true(all(res %% 1 ==0))
})

test_that("tsourakakis_2014_approximate gives the correct solution", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2014_approximate(g)
  
  expect_true(all(res == (n/2+1):n))
})

test_that("tsourakakis_2014_approximate works for detecting the larger clique", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2014_approximate(g)
  
  expect_true(all(res == as.character((n/3+1):n)))
})

test_that("tsourakakis_2014_approximate works a random graph", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2014_approximate(g)
  
  expect_true(is.numeric(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
  expect_true(length(res) > 0)
})

test_that("tsourakakis_2014_approximate works with a core set", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n-10, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  edge_mat <- cbind(edge_mat, rbind(20:29, 21:30))
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2014_approximate(g, core_set = c(10:14, n))
  
  expect_true(n %in% res)
})

#################3

## anderson_2009 is correct
test_that("anderson_2009 works", {
  n <- 10
  g <- igraph::make_ring(n)
  res <- anderson_2009(g)
  
  expect_true(is.numeric(res))
  expect_true(length(res) <= n)
  expect_true(length(unique(res)) == length(res))
  expect_true(all(res > 0))
  expect_true(all(res <= n))
  expect_true(all(res %% 1 ==0))
})

test_that("anderson_2009 can find large clique", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- anderson_2009(g)
  
  expect_true(all(res == (n/2+1):n))
})

test_that("anderson_2009 works with an empty graph", {
  n <- 20
  g <- igraph::graph.empty(n = n, directed = F)
  res <- anderson_2009(g)
  
  expect_true(length(res) == 0)
})

test_that("anderson_2009 works for detecting the larger clique", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- anderson_2009(g)
  
  expect_true(all(res == as.character((n/3+1):n)))
})

test_that("anderson_2009 works a random graph", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- anderson_2009(g)
  
  expect_true(is.numeric(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
  expect_true(length(res) > 0)
})

test_that("anderson_2009 works with a core set", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n-10, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  edge_mat <- cbind(edge_mat, rbind(20:29, 21:30))
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- anderson_2009(g, core_set = c(10:14, n))
  
  expect_true(n %in% res)
})

#######################

## .chen_check_density is correct
test_that(".chen_check_density works", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  node_set <- as.character((n/2+1):n)
  
  res <- .chen_check_density(g, node_set)
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".chen_check_density gets the right value", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  node_set <- as.character((n/2+1):n)
  
  res <- .chen_check_density(g, node_set)
  
  expect_true(abs(res - floor(0.975*ncol(tmp))/ncol(tmp)) <= 1e-6)
})

## .form_c_matrix is correct
test_that(".form_c_matrix works", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  node_set <- as.character((n/2+1):n)
  
  res <- .form_c_matrix(g)
  
  expect_true(is.matrix(res))
  expect_true(ncol(res) == 3)
  expect_true(nrow(res) == choose(n,2))
})

## .chen_separate is correct
test_that(".chen_separate works", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)
  edge_mat <- cbind(c(1,n), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))], tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))] + n/2)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  c_mat <- .form_c_matrix(g)
  res <- .chen_separate(g, c_mat, check = T)
  
  expect_true(length(res) == 4)
  expect_true(class(res$g1) == "igraph")
  expect_true(class(res$g2) == "igraph")
  expect_true(is.matrix(res$c_matrix1))
  expect_true(is.matrix(res$c_matrix2))
  expect_true(all(dim(res$c_matrix1) == c(choose(igraph::vcount(res$g1),2), 3)))
  expect_true(all(dim(res$c_matrix2) == c(choose(igraph::vcount(res$g2),2), 3)))
})

test_that(".chen_separate is correct", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)
  edge_mat <- cbind(c(1,n), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))], tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))] + n/2)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  c_mat <- .form_c_matrix(g)
  res <- .chen_separate(g, c_mat, check = T)
  
  expect_true(all(igraph::V(res$g1)$name == as.character(1:(n/2))) | all(igraph::V(res$g1)$name == as.character((n/2+1):n)) )
  expect_true(all(igraph::V(res$g2)$name == as.character(1:(n/2))) | all(igraph::V(res$g2)$name == as.character((n/2+1):n)) )
  expect_true(all(res$c_matrix1[,1] %in% igraph::V(res$g1)$name))
  expect_true(all(res$c_matrix1[,2] %in% igraph::V(res$g1)$name))
  expect_true(all(res$c_matrix2[,1] %in% igraph::V(res$g2)$name))
  expect_true(all(res$c_matrix2[,2] %in% igraph::V(res$g2)$name))
})

## chen_2010 is correct
test_that("chen_2010 works", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)
  edge_mat <- cbind(c(1,n), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))], tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))] + n/2)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- chen_2010(g, threshold = 0.95)
  
  expect_true(all(res == as.character(1:(n/2))) | all(res == as.character((n/2+1):n)))
})

test_that("chen_2010 works for high node disparity", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp[,sample(ncol(tmp), size = floor(0.975*ncol(tmp)))])
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- chen_2010(g, threshold = 0.95)
  
  expect_true(all(res == as.character((n/2+1):n)))
})

test_that("chen_2010 works for detecting the larger clique", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- chen_2010(g, threshold = 0.95)
  
  expect_true(all(res == as.character((n/3+1):n)))
})

test_that("chen_2010 works a random graph", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- chen_2010(g)
  
  expect_true(is.numeric(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
  expect_true(length(res) > 0)
})

test_that("chen_2010 works with a core set", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n-10, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  edge_mat <- cbind(edge_mat, rbind(20:29, 21:30))
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- chen_2010(g, core_set = c(10:14, n))
  
  expect_true(n %in% res)
})


###############

## .tsourakakis_initialize is correct
test_that(".tsourakakis_initialize works", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  res <- .tsourakakis_initialize(g)
  
  expect_true(is.character(res))
  expect_true(length(res) == 1)
})

## .tsourakakis_obj is correct
test_that(".tsourakakis_obj works", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  res <- .tsourakakis_obj(g, 0.95, as.character(11:20))
  
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".tsourakakis_obj can be negative", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  res <- .tsourakakis_obj(g, 0.95, as.character(1:10))
  
  expect_true(res < 0)
  expect_true(res == 9 - 0.95*10*9/2)
})

test_that(".tsourakakis_obj can be positive", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  res <- .tsourakakis_obj(g, 0.95, as.character(11:20))
  
  expect_true(res > 0)
  expect_true(abs(res -10*9/2 * (1-0.95)) <= 1e-6)
})

## .find_candidate is correct
test_that(".find_candidate can add", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  
  res <- .find_candidate(g, 0.95, node_set = as.character(16:20), node_candidate = as.character(1:15),
                         den_org = 0.5)
  
  expect_true(length(res) == 6)
  expect_true(length(intersect(res, as.character(16:20))) == 5)
})

test_that(".find_candidate can remove", {
  set.seed(10)
  n <- 20
  tmp <- utils::combn(n/2, 2)+n/2
  edge_mat <- cbind(rbind(1:(n/2-1), 2:(n/2)), tmp)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  
  den_target <- .tsourakakis_obj(g, 0.95, as.character(16:20))
  
  res <- .find_candidate(g, 0.95, node_set = as.character(c(1,16:20)), node_candidate = as.character(c(1,16:20)),
                         den_org = 0, remove = T)
  
  expect_true(length(res) == 5)
  expect_true(length(intersect(res, as.character(c(1,16:20)))) == 5)
})

## tsourakakis_2013 is correct
test_that("tsourakakis_2013 works", {
  set.seed(10)
  n <- 30
  tmp1 <- utils::combn(n/3, 2)
  tmp2 <- utils::combn(2*n/3, 2)
  edge_mat <- cbind(c(1,n), tmp1[,sample(ncol(tmp1), size = floor(0.975*ncol(tmp1)))], tmp2[,sample(ncol(tmp2), size = floor(0.975*ncol(tmp2)))] + n/3)
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  igraph::V(g)$name <- 1:n
  res <- tsourakakis_2013(g)
  
  expect_true(is.numeric(res))
  expect_true(length(res) > n/3)
})

test_that("tsourakakis_2013 works a random graph", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2013(g)
  
  expect_true(is.numeric(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
  expect_true(length(res) > 0)
})

test_that("tsourakakis_2013 works with a core set", {
  set.seed(10)
  n <- 30
  combn_mat <- utils::combn(n-10, 2)
  edge_mat <- combn_mat[,sample(1:ncol(combn_mat), round(ncol(combn_mat)/2))]
  edge_mat <- cbind(edge_mat, rbind(20:29, 21:30))
  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  res <- tsourakakis_2013(g, core_set = c(10:14, n))
  
  expect_true(n %in% res)
})


