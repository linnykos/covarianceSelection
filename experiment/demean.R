set.seed(10); n <- 10
x <- cbind(rnorm(n, mean = 3), rnorm(n, mean = 3))
x1 <- scale(x, center = T, scale = F)
one <- rep(1, n); decomp <- base::qr(one)
Q <- base::qr.Q(decomp, complete=TRUE); Q <- Q[,2:ncol(Q),drop = F]
x2 <- t(Q)%*%x

xlim <- c(min(x[,1], x1[,1], x2[,1]), max(x[,1], x1[,1], x2[,1]))
ylim <- c(min(x[,2], x1[,2], x2[,2]), max(x[,2], x1[,2], x2[,2]))

par(mar = rep(1,4)); plot(x, xlim = xlim, ylim = ylim, pch = 16)
points(x1, pch = 16, col = "red"); points(x2, pch = 16, col = "green")

x1[-1,] - x2
apply(x2, 2, mean) - apply(x1[-1,], 2, mean)
apply(x[-1,],2,sd); apply(x2,2,sd)

t(Q)[1,]%*%x-colMeans(x2)
x[2,]-colMeans(x)-colMeans(x1[-1,])

