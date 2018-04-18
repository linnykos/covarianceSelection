set.seed(10)
n <- 1000
dat <- stats::rnorm(n)
#dat <- c(stats::rnorm(n/2), stats::rnorm(n/2, mean = 5))

proportion <- seq(0.1,1,length.out = 50)
mean_val <- mean(dat)
trials <- 1000

vec <- rep(NA, length(proportion))
for(i in 1:length(proportion)){
  tmp <- sapply(1:trials, function(x){
    subsamp <- dat[sample(floor(proportion[i]*n))]
    ci_ss <- t.test(subsamp, conf.level = 0.99)$conf.int
    ifelse(ci_ss[1] <= mean_val & ci_ss[2] >= mean_val, TRUE, FALSE)
  })
  vec[i] <- sum(tmp)
}

vec/trials


