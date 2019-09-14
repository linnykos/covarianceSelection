rm(list=ls())
set.seed(10)
trials <- 1000
d <- 3; n <- 50
i <- 1
set.seed(10*i)
dat_list <- lapply(1:5, function(x){
  mat <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
})
alpha = 0.05
return_pvalue = T
cores = 1
verbose = T

doMC::registerDoMC(cores = cores)

dat_list <- lapply(dat_list, scale, center = T, scale = F)
len <- length(dat_list)
combn_mat <- utils::combn(len, 2)
idx_all <- rep(TRUE, ncol(combn_mat))

diag_idx <- which(lower.tri(diag(ncol(dat_list[[1]])), diag = T))
num_list <- lapply(dat_list, function(x){.compute_sigma(x, diag_idx)})
denom_list <- .compute_all_denom(dat_list, num_list, diag_idx)

t_vec <- .compute_all_test_stat(num_list, denom_list, combn_mat = combn_mat, squared = T)

if(verbose)  print(paste0("Starting to run heavy parallel computation: ", Sys.time()))

func <- function(i){
  if(verbose && i %% floor(trials/10) == 0) cat('*')
  set.seed(round*10*i)
  noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
  
  remaining_pairs <- which(idx_all)
  combn_short <- combn_mat[, remaining_pairs, drop = F]
  if(any(is.na(combn_short))) return(Inf)
  
  remaining_idx <- unique(as.vector(combn_short))
  boot_num_list <- .compute_all_numerator_bootstrap(dat_list, noise_list, num_list, diag_idx,
                                               remaining_idx = 1:length(dat_list))
  
  boot_t_vec <- .compute_all_test_stat(boot_num_list, denom_list, combn_mat = combn_mat)
  
  if(round == 1 & return_pvalue){
    list(val = max(abs(boot_t_vec)), boot_t_vec = boot_t_vec)
  } else {
    list(val = max(abs(boot_t_vec)), boot_t_vec = NA)
  }
}

round <- 1
round_1_boot_t_vec <- NA

if(sum(idx_all) == 0) break()

i <- 0 #debugging purposes
res <- foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i))
t_boot <- sapply(res, function(x){x$val})
if(round == 1 & return_pvalue) round_1_boot_t_vec <- sapply(res, function(x){x$boot_t_vec})

cutoff <- stats::quantile(t_boot, 1-alpha)
idx <- intersect(which(abs(t_vec) >= cutoff), which(idx_all))

idx_all[idx] <- FALSE

stopifnot(length(t_vec) == nrow(round_1_boot_t_vec))
pval <- sapply(1:length(t_vec), function(i){
  length(which(abs(round_1_boot_t_vec[i,]) > abs(t_vec[i])))/ncol(round_1_boot_t_vec)
})

############################3

# extensive side-by-side comparison
x <- dat_list[[1]]
y <- dat_list[[2]]
diag_idx <- which(lower.tri(diag(ncol(x)), diag = T))

n1 <- nrow(x); n2 <- nrow(y)

num_x <- .compute_sigma(x, diag_idx); num_y <- .compute_sigma(y, diag_idx)
denom_x <- .compute_variance(x, num_x, diag_idx); denom_y <- .compute_variance(y, num_y, diag_idx)
t_org <- .compute_covStat(num_x, num_y, denom_x, denom_y)

func <- function(i){
  if(verbose && i %% floor(trials/10) == 0) cat('*')
  set.seed(round*10*i)
  noise_list <- lapply(dat_list, function(x){stats::rnorm(nrow(x))})
  
  g_x <- noise_list[[1]]; g_y <- noise_list[[2]]
  boot_num_x <- .compute_bootSigma(x, g_x, num_x, diag_idx)
  boot_num_y <- .compute_bootSigma(y, g_y, num_y, diag_idx)
  .compute_covStat(boot_num_x, boot_num_y, denom_x, denom_y)
}

if(is.na(cores)){
  t_boot <- numeric(trials)
  for(i in 1:trials){ t_boot[i] <- func(i) }
} else {
  i <- 0 #debugging purposes
  t_boot <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i)))
}

length(which(abs(t_boot) >= abs(t_org)))/trials
