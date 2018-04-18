pthres <- 0.05
lambda <- 0.12 # from https://arxiv.org/pdf/1506.00728.pdf, top of page 23
num_target <- 246 # from https://arxiv.org/pdf/1506.00728.pdf, bottom of page 23

if(verbose) print("Forming graph")
prec_mat <- covarianceSelection::graphicalModel_form(dat, filename_func, lambda = lambda,
                                                num_primary = num_primary,
                                                cores = cores, verbose = verbose)
if(verbose) print(paste0(Sys.time(), ": Finished forming graph"))
sigma_mat <- stats::cov(dat)
d <- ncol(dat)
tmp <- prec_mat; tmp[lower.tri(tmp, diag = T)] <- 0
edges <- which(abs(tmp) > 1e-6, arr.ind = T)

if(verbose) print(paste0(Sys.time(), ": Edge testing. Graph has ", nrow(prec_mat), " nodes and ",
                         nrow(edges), " edges."))
adj_gene <- matrix(0, d, d)
for(i in 1:nrow(edges)){
  adj_gene[edges[i,1], edges[i,2]] <- 1; adj_gene[edges[i,2], edges[i,1]] <- 1
}
colnames(adj_gene) <- tada$Gene; rownames(adj_gene) <- tada$Gene

seedindex <- rep(0, ncol(dat))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_res <- covarianceSelection::hmrf(tada$pval.TADA, adj_gene, seedindex, pthres = pthres)
report <- covarianceSelection::report_results(tada$Gene, 1-hmrf_res$post, tada$pval.TADA, hmrf_res$Iupdate)
cutoff <- sort(report$FDR, decreasing = FALSE)[num_target]
autism_genes <- report$Gene[which(report$FDR <= cutoff)]

rm(list = c("sigma_mat", "d", "tmp", "edges", "seedindex"))

save.image(file = paste0(save_filepath, "step4_res", additional_name, ".RData"))
