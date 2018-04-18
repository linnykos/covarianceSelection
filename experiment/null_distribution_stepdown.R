trials <- 5
n <- 500; d <- 5; k <- 2; p <- 3
cov_list <- list(diag(d), huge::huge.generator(n = 2, d = d, graph = "random", verbose = F)$sigma)

generate_data <- function(seed, n, d, k, p){
  set.seed(seed)

  dat_list <- vector("list", p*k)
  for(i in 1:k){
    for(j in 1:p){
      dat_list[[(i-1)*p+j]] <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = cov_list[[i]])
    }
  }

  dat_list
}

combn_mat <- combn(p*k, 2)
bool_vec <- apply(combn_mat, 2, function(x){
  if(floor((x[1]-.5)/p) == floor((x[2]-.5)/p)) return(TRUE) else return(FALSE)
})
correct_idx <- which(bool_vec)

func <- function(trial){
  dat_list <- generate_data(trial, n = n, d = d, k = k, p = p)
  set.seed(trial)
  longitudinalGM::stepdown(dat_list, trials = 1000, cores = 1, alpha = 1)
}

doMC::registerDoMC(cores = 3)
res <- foreach::"%dopar%"(foreach::foreach(trial = 1:trials), func(trial))
save.image("../experiment/null_distribution_stepdown.RData")

load("../experiment/null_distribution_stepdown.RData")
count_vec <- sapply(res, function(x){length(which(x %in% correct_idx))})
length(which(count_vec == length(correct_idx)))
0.95*trials - 3*sqrt(0.95*0.05*trials)
