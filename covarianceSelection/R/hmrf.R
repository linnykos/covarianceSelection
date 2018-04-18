#partial likelihood function
.partial_likelihood <- function(b, c, graph_term, i_vec){
  d <- length(i_vec)

  fv.fun <- function(i){
    new1 <- exp(b*i_vec[i]+c*i_vec[i]*graph_term[i])
    new2 <- exp(b*(1-i_vec[i])+c*(1-i_vec[i])*graph_term[i])
    fvalue <- log(new1/(new1+new2))
  }

  sum(sapply(1:d,fv.fun))
}

##estimate Ising model parameters b,c
.optimize_bc <- function(func, adj, i_vec, times, tol = 1e-5){
  b <- 0; c <- 1

  graph_term <- i_vec%*%adj
  for (k in 1:times){
    b_new <- stats::optimize(func, c(-2,2), c = c, graph_term = graph_term,
                             i_vec = i_vec, maximum = T)$maximum
    c_new <- stats::optimize(func, c(0.2,2), b = b_new, graph_term = graph_term,
                             i_vec = i_vec, maximum = T)$maximum
    if (abs(b_new-b) < tol & abs(c_new-c) < tol){
      break()
    }

    b <- b_new; c <- c_new
  }

  list(b = b, c = c)
}


#' Hidden markov random field
#'
#' Use a Ising model to update the pvalues according to the graphical structure.
#' There are two groups, indices in group 0 have mean 0. Indices in group 1 have
#' mean greater than 0.
#'
#' @param pval a vector of d p-values, between 0 and 1.
#' @param adj a dxd adjancey matrix (0,1 entries)
#' @param seedindex a (0,1) vector of length d, where 1 means the gene is in group 1.
#' @param pthres threshold for p-values to serve as an initialization
#' @param iter number of iterations
#' @param verbose boolean
#' @param tol tolerance for values that should be treated as zero in absolute value
#'
#' @return a list
#' @export
hmrf <- function(pval, adj, seedindex, pthres = 0.05, iter = 100,
                 verbose = FALSE, tol = 1e-3){
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

  for (iteri in 1:iter){

    if(verbose) print(paste("Start iteration: ", iteri))
    res <- .optimize_bc(.partial_likelihood, adj, i_vec, 20)
    b_new <- res$b; c_new <- res$c

    if (abs(c-c_new)<tol & abs(b-b_new)<tol)  break()
    b <- b_new; c <- c_new

    for (i in 1:d){
      i_vec_tmp <- i_vec; i_vec_tmp[i] <- 1-i_vec[i]

      new1 <- exp(b*i_vec[i] + c*i_vec[i]*adj[i,]%*%i_vec)
      new2 <- exp(b*(1-i_vec[i]) + c*(1-i_vec[i])*adj[i,]%*%i_vec_tmp)
      p1 <- stats::dnorm(z[i], mu1*i_vec[i], sqrt(sigmas)) * new1/(new1+new2)
      p2 <- stats::dnorm(z[i], mu1*(1-i_vec[i]), sqrt(sigmas)) * new2/(new1+new2)

      if (i_vec[i] == 1){
        posterior[i] <- p1/(p1+p2)
      } else {
        posterior[i] <- p2/(p1+p2)
      }

      if (p2 > p1){ i_vec[i] <- 1-i_vec[i] }
      if (seedindex[i] != 0){ i_vec[i] <- 1 }
    }

    mu1 <- sum(posterior[seedindex==0]*z[seedindex==0])/sum(posterior[seedindex==0])
    sigmas2 <- sum(posterior[seedindex==0]*(z[seedindex==0]-mu1)^2)/sum(posterior[seedindex==0])
    sigmas1 <- sum((1-posterior[seedindex==0])*(z[seedindex==0])^2)/sum(1-posterior[seedindex==0])
    sigmas <- (sigmas1*sum(posterior[seedindex==0])+sigmas2*sum(1-posterior[seedindex==0])) / length(posterior)

    if(verbose) print(paste0("Iteration: ",iteri," has ", sum(i_vec)," genes set with Iupdate = 1."))
  }

  #un-permute the elements
  i_vec <- i_vec[order(idx)]; posterior <- posterior[order(idx)]

  list(Iupdate = i_vec, post = posterior, b = b, c = c, mu1 = mu1, sigmas = sigmas)
}
