rm(list=ls())

ncores <- 20
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
devtools::install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)
library(org.Hs.eg.db)

verbose <- T
save_filepath <- "/raid6/Kevin/covarianceSelection/results"
tmp <- 1-1e-4/2
filepath_suffix <- paste0("_", tmp)

source("../main/step0_loading.R")
source("../main/step1_screening.R")
source("../main/step2_nodawn_analysis.R")
source("../main/step3_pfc35_analysis.R")
source("../main/step4_alldata_analysis.R")
source("../main/step5_subjectselection.R")
source("../main/step6_our_analysis.R")
source("../main/step7_results.R")
