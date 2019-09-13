d <- 5
n <- 20
set.seed(10)
x <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
y <- MASS::mvrnorm(n = n, mu = rep(0,d), Sigma = diag(d))
# res <- cai_test(x,y, trials = 50)

trials = 100
cores = 1

if(ncol(x) != ncol(y)) stop("x and y have different number of dimensions")
if(!is.matrix(x) | !is.numeric(x)) stop("x is not a numeric matrix")
if(!is.matrix(y) | !is.numeric(y)) stop("y is not a numeric matrix")

doMC::registerDoMC(cores = cores)
diag_idx <- which(lower.tri(diag(ncol(x)), diag = T))

n1 <- nrow(x); n2 <- nrow(y)

num_x <- .compute_sigma(x, diag_idx); num_y <- .compute_sigma(y, diag_idx)
denom_x <- .compute_variance(x, num_x, diag_idx); denom_y <- .compute_variance(y, num_x, diag_idx)
t_org <- .compute_covStat(num_x, num_y, denom_x, denom_y)

func <- function(i){
  set.seed(i*10)
  g_x <- stats::rnorm(n1); g_y <- stats::rnorm(n2)
  boot_num_x <- .compute_bootSigma(x, g_x, num_x, diag_idx)
  boot_num_y <- .compute_bootSigma(y, g_y, num_y, diag_idx)
  .compute_covStat(boot_num_x, boot_num_y, denom_x, denom_y)
}
