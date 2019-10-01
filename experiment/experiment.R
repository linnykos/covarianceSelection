rm(list=ls())
ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)
load("/raid6/Kevin/covarianceSelection/results/step2_nodawn_analysis.RData")

#####
set.seed(10)
if(verbose) print(paste0(Sys.time(), "Start of step 2: Naive analysis"))

selected_idx <- grep("PFC\\.[3-5]", names(dat_list))
dat_pfc35 <- do.call(rbind, dat_list[selected_idx]) # 107 x 3438
dat_pfc35 <- scale(dat_pfc35)

# estimate graphical model on PFC35 using cross-validated lasso for neighborhood selection
res <- covarianceSelection:::graphicalModel_range(dat_pfc35, lambda_min = 0.01, lambda_max = 0.35, verbose = T) 
save.image("../experiment/tmp.RData")