cores <- 10
set.seed(10)

library(devtools)
#install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)
library(foreach)
library(doMC)
library(hash)
library(glmnet)
library(fdrtool)

verbose = F
save_filepath = "../results/"
