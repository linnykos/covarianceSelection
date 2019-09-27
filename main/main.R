rm(list=ls())

ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
devtools::install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)

verbose <- T
save_filepath <- "/raid6/Kevin/covarianceSelection/results"

source("../main/step0_loading.R")
source("../main/step1_screening.R")
source("../main/step2_nodawn_analysis.R")
source("../main/step3_pfc35_analysis.R")
source("../main/step4_alldata_analysis.R")
source("../main/step5_subjectselection.R")


# load("../results/step2_naive_analysis.RData")
# 
# save_filepath <- "/raid6/Kevin/covarianceSelection/results"
# ncores <- 20
# source("../main/step3_subjectselection.R")