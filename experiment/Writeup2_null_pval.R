rm(list=ls())
trials <- 300
d <- 3; n <- 50

p_null <- numeric(trials)
for(i in 1:trials){
  set.seed(10*i)
  x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
  p_null[i] <- cai_test(x,y, trials = 300)
}

par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")

##################

set.seed(10)
trials <- 300
d <- 3; n <- 500

res <- lapply(1:trials, function(i){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(10*i)
  dat_list <- lapply(1:5, function(x){MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))})
  stepdown(dat_list, trials = 1000, return_pvalue = T)
})

p_null <- unlist(lapply(res, function(x){x$pval}))
par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")

# what if you plot things marginally
p_null <- unlist(lapply(res, function(x){x$pval[1]}))
par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")

#########################

set.seed(10)
trials <- 300
d <- 3; n <- 50

res <- lapply(1:trials, function(i){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(10*i)
  dat_list <- lapply(1:5, function(x){MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))})
  cai_test(dat_list[[1]], dat_list[[2]], trials = 1000)
})

p_null <- unlist(res)
par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")

######################

set.seed(10)
trials <- 50
d <- 3; n <- 1000

res <- lapply(1:trials, function(i){
  print(i)
  set.seed(10*i)
  dat_list <- lapply(1:5, function(x){MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))})
  obj <- stepdown_path(dat_list, trials = 1000, iterations = 10)
  stepdown_choose(obj, alpha = 0.05, return_pvalue = T)
})

p_null <- unlist(lapply(res, function(x){x$pval}))
par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")


# what if you plot things marginally
p_null <- unlist(lapply(res, function(x){x$pval[1]}))
par(mfrow = c(1,2))
plot(sort(p_null), seq(0,1,length.out=length(p_null)))
lines(c(0,1), c(0,1), col = "red")
hist(p_null, breaks = 10, col = "gray")
