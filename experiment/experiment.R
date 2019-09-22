rm(list=ls())
load("../results/step3_alldata_analysis.RData")

set.seed(10)
seedindex <- rep(0, ncol(adj_all))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_all <- covarianceSelection::hmrf(tada$pval.TADA, adj_all, seedindex, pthres = pthres)
####################

pval = tada$pval.TADA
adj = adj_all
iter = 100
verbose = FALSE
tol = 1e-3
pthres = 0.05
 
stopifnot(length(pval) == length(seedindex), all(dim(adj) == length(pval)))

#permute all the entries to avoid emphasis on current order
d <- length(pval)
idx <- sample(1:d)
pval <- pval[idx]; adj <- adj[idx,idx]; seedindex <- seedindex[idx]

z <- stats::qnorm(1-pval)
i_vec <- as.numeric(pval<pthres)
b <- -Inf; c <- Inf

mu1 <- mean(z[i_vec==1 & seedindex==0])
sigmas <- stats::var(z)
posterior <- rep(0,d)

# res <- .optimize_bc(.partial_likelihood, adj, i_vec, 20)
##############

func <- .partial_likelihood
tol = 1e-5
times <- 20
b <- 0; c <- 0

graph_term <- i_vec%*%adj
for (k in 1:times){
  b_new <- stats::optimize(func, c(-20, 2), c = c, graph_term = graph_term,
                           i_vec = i_vec, maximum = T)$maximum
  c_new <- stats::optimize(func, c(-2, 2), b = b_new, graph_term = graph_term,
                           i_vec = i_vec, maximum = T)$maximum
  if (abs(b_new-b) < tol & abs(c_new-c) < tol){
    break()
  }
  print(paste0("b: ", b_new, "// c: ", c_new))
  
  b <- b_new; c <- c_new
}
