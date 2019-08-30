pthres <- 0.05
num_target <- 200

#####

if(verbose) print(paste0(Sys.time(), "Start of step 2: Naive analysis"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3438
dat_pfc35 <- scale(dat_pfc35)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
prec_mat_naive <- covarianceSelection::graphicalModel(dat_pfc35, lambda = "lambda.min", verbose = T) 

# grab the edges
d <- ncol(dat_pfc35)
tmp <- prec_mat_naive; tmp[lower.tri(tmp, diag = T)] <- 0
edges_naive <- which(abs(tmp) > 1e-6, arr.ind = T)

# form adjacency matrix
if(verbose) print(paste0(Sys.time(), ": Graph has ", nrow(prec_mat_naive), " nodes and ",
                         nrow(edges_naive), " edges."))
adj_gene <- matrix(0, d, d)
for(i in 1:nrow(edges_naive)){
  adj_gene[edges_naive[i,1], edges_naive[i,2]] <- 1; adj_gene[edges_naive[i,2], edges_naive[i,1]] <- 1
}
colnames(adj_gene) <- tada$Gene; rownames(adj_gene) <- tada$Gene

# run the HMRF
set.seed(10)
seedindex <- rep(0, ncol(dat_pfc35))
seedindex[which(tada$dn.LoF >= 3)] <- 1

if(verbose) print(paste0(Sys.time(), ": HMRF"))
hmrf_res <- covarianceSelection::hmrf(tada$pval.TADA, adj_gene, seedindex, pthres = pthres)
report <- covarianceSelection::report_results(tada$Gene, 1-hmrf_res$post, tada$pval.TADA, hmrf_res$Iupdate)
cutoff <- sort(report$FDR, decreasing = FALSE)[num_target]
autism_genes <- report$Gene[which(report$FDR <= cutoff)]


save.image(file = paste0(save_filepath, "/step2_naive_analysis.RData"))