n <- 5; d <- 5
one <- matrix(1, ncol = 1, nrow = n)

decomp <- base::qr(one)
Q <- base::qr.Q(decomp, complete=TRUE)

set.seed(10)
mat <- matrix(rnorm(n*d), n, d)
z <- Q%*%mat

mat[2,]%*%mat[2,] - z[2,]%*%z[2,]
mat[3,]%*%mat[3,] - z[3,]%*%z[3,]
mat[2,]%*%mat[3,] - z[2,]%*%z[3,]

l2norm <- function(x){sqrt(sum(x^2))}
l2norm(mat[2,]-mat[3,])
l2norm(z[2,]-z[3,])

mat %*% t(mat)
mat[2,]%*%mat[2,]
Q[2,]%*%t(Q[2,])

trace <- function(mat){sum(diag(mat))}
P = Q[2,]%*%t(Q[2,])
D = mat%*%t(mat)
trace(P%*%D)
z[2,]%*%z[2,]

Q[2,]-Q[3,]

##################
#reparameterizing Q
project_out <- function(vec1, vec2){
  d <- length(vec1)
  vec2 <- (diag(d) - vec1%*%t(vec1))%*%vec2
  vec2/l2norm(vec2)
}

project_out_many <- function(list_vec, vec2){
  k <- length(list_vec)
  for(i in 1:k){
    vec2 <- project_out(list_vec[[i]], vec2)
  }

  vec2
}

reparameterize_basis <- function(mat){
  d <- nrow(mat)
  one <- rep(1, d)
  list_vec <- vector("list", 1); list_vec[[1]] <- one/l2norm(one)
  for(i in 1:(d-1)){
    vec2 <- project_out_many(list_vec, mat[,i])
    list_vec[[i+1]] <- vec2
  }

  do.call(cbind, list_vec)
}

check_valid_q <- function(mat, tol = 1e-6){
  d <- nrow(mat)
  one <- matrix(1, ncol = 1, nrow = d)
  r <- c(l2norm(one), rep(0, d-1))
  tmp <- mat %*% r

  stopifnot(sum(abs(tmp - one)) <= tol)
  stopifnot(sum(abs(mat%*%t(mat) - diag(d))) <= tol)
  stopifnot(sum(abs(t(mat)%*%mat - diag(d))) <= tol)

  TRUE
}

set.seed(40)
mat <- matrix(rnorm(100),10,10)
mat <- eigen(mat%*%t(mat))$vectors
mat <- reparameterize_basis(mat)
#mat%*%t(mat) #for checking
Q <- mat
check_valid_q(Q)

set.seed(10); n <- 10
x <- cbind(rnorm(n, mean = 3), rnorm(n, mean = 3))
x1 <- scale(x, center = T, scale = F)
Q <- Q[,2:ncol(Q),drop = F]
x2 <- t(Q)%*%x

xlim <- c(min(x[,1], x1[,1], x2[,1]), max(x[,1], x1[,1], x2[,1]))
ylim <- c(min(x[,2], x1[,2], x2[,2]), max(x[,2], x1[,2], x2[,2]))

par(mar = rep(1,4)); plot(x, xlim = xlim, ylim = ylim, pch = 16, asp = T)
points(x1, pch = 16, col = "red"); points(x2, pch = 16, col = "green")

