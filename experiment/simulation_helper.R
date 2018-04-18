clean_matrix <- function(mat){
  diag(mat) <- 1
  mat <- (mat + t(mat))/2
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0
  mat <- eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
  diag(1/diag(mat)^(1/2))%*%mat%*%diag(1/diag(mat)^(1/2))
}

#fixed design with four clusters
generate_matrix <- function(d = 60, strength = 0){
  stopifnot(d %% 4 == 0)

  d4 <- d/4
  vec1 <- c(.6 + strength*.4, strength, strength, .7 + strength*.3)
  vec2 <- c(strength, .9 + strength*.1, .5 + strength*.5, strength)
  vec3 <- c(strength, .5 + strength*.5, .4 + strength*.6, .7 + strength*.3)
  vec4 <- c(.7 + strength*.3, strength, .7 + strength*.3, .9 + strength*.1)

  cov_mat <- lapply(1:4, function(x){
    if(x == 1){
      matrix(rep(rep(vec1, each = d4), times = d4), ncol = d4, nrow = d)
    } else if (x == 2){
      matrix(rep(rep(vec2, each = d4), times = d4), ncol = d4, nrow = d)
    } else if (x == 3){
      matrix(rep(rep(vec3, each = d4), times = d4), ncol = d4, nrow = d)
    } else {
      matrix(rep(rep(vec4, each = d4), times = d4), ncol = d4, nrow = d)
    }
  })

  cov_mat <- do.call(cbind, cov_mat)
}

error_function <- function(cluster1, cluster2){
  stopifnot(length(cluster1) == length(cluster2))
  d <- length(cluster1)

  combn_mat <- combn(d, 2)
  m <- ncol(combn_mat)
  res <- apply(combn_mat, 2, function(x){
    bool1 <- cluster1[x[1]] == cluster1[x[2]]
    bool2 <- cluster2[x[1]] == cluster2[x[2]]
    if(bool1 == bool2) TRUE else FALSE
  })
  sum(res)/m
}
