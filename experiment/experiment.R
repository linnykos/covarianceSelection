set.seed(10)
d <- 12
pval <- stats::runif(d)
adj <- as.matrix(huge::huge.generator(n = 50, d = d, graph = "scale-free", verbose = F)$theta)
seedindex <- rep(0, d)
seedindex[order(pval, decreasing = F)[1:2]] <- 1

pthres = 0.5
iter = 2
verbose = FALSE
tol = 1e-3

##

stopifnot(length(pval) == length(seedindex), all(dim(adj) == length(pval)))

#permute all the entries to avoid emphasis on current order
d <- length(pval)
idx <- sample(1:d)
pval <- pval[idx]; adj <- adj[idx,idx]; seedindex <- seedindex[idx]

z <- stats::qnorm(1-pval)
i_vec <- as.numeric(pval<pthres)
b <- Inf; c <- Inf

non_seedidx <- (seedindex==0)
mu1 <- mean(z[i_vec==1 & non_seedidx])
sigmas1 <- (stats::sd(z[i_vec == 0]))^2
sigmas2 <- (stats::sd(z[i_vec == 1 & non_seedidx]))^2
posterior <- rep(0,d)

for (iteri in 1:iter){
  
  if(verbose) print(paste("Start iteration: ", iteri))
  res <- .optimize_bc(.partial_likelihood, adj, i_vec, 20)
  b_new <- res$b; c_new <- res$c
  
  if (abs(c-c_new)<tol & abs(b-b_new)<tol)  break()
  b <- b_new; c <- c_new
  
  for (i in 1:d){
    i_vec_tmp <- i_vec; i_vec_tmp[i] <- 1-i_vec[i]
    
    new1 <- b*i_vec[i] + c*i_vec[i]*adj[i,]%*%i_vec
    new2 <- b*(1-i_vec[i]) + c*(1-i_vec[i])*adj[i,]%*%i_vec_tmp
    p1 <- stats::dnorm(z[i], mu1*i_vec[i], sqrt(sigmas2 * i_vec[i] + sigmas1 * (1 - i_vec[i])))/(1 + exp(new2-new1))
    p2 <- stats::dnorm(z[i], mu1*(1-i_vec[i]), sqrt(sigmas2 * (1 - i_vec[i]) + sigmas1 * i_vec[i]))/(1 + exp(new1-new2))
    
    if (i_vec[i] == 1){
      posterior[i] <- p1/(p1+p2)
    } else {
      posterior[i] <- p2/(p1+p2)
    }
    
    if (p2 > p1){ i_vec[i] <- 1-i_vec[i] }
    if (seedindex[i] != 0){ i_vec[i] <- 1 }
  }
  
  mu1 <- sum(posterior[non_seedidx]*z[non_seedidx])/sum(posterior[non_seedidx])
  sigmas2 <- sum(posterior[non_seedidx]*(z[non_seedidx]-mu1)^2)/sum(posterior[non_seedidx])
  sigmas1 <- sum((1-posterior[non_seedidx])*(z[non_seedidx])^2)/sum(1-posterior[non_seedidx])
  sigmas <- (sigmas1*sum(posterior[non_seedidx]) + sigmas2*sum(1-posterior[non_seedidx])) / length(posterior[non_seedidx])
  sigmas1 <- sigmas; sigmas2 <- sigmas
  
  if(verbose) print(paste0("Iteration: ",iteri," has ", sum(i_vec)," genes set with Iupdate = 1."))
}