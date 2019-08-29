rm(list=ls())

cores <- 10
set.seed(10)

library(devtools)
#install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)

verbose <- F
save_filepath <- "../results/"

source("../main/step0_loading.R")