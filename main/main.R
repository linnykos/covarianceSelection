rm(list=ls())

ncores <- 15
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
#install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)

verbose <- T
save_filepath <- "../results/"

source("../main/step0_loading.R")
source("../main/step1_screening.R")
source("../main/step2_naive_analysis.R")