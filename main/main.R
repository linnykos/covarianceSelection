rm(list=ls())

ncores <- 10
set.seed(10)
doMC::registerDoMC(cores = ncores)

library(devtools)
devtools::install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection", force = T)
library(covarianceSelection)
library(org.Hs.eg.db)

verbose <- T
save_filepath <- "/raid6/Kevin/covarianceSelection/results"
filepath_suffix <- ""

source("../main/step0_loading.R")
source("../main/step1_screening.R")
source("../main/step2_nodawn_analysis.R")
source("../main/step3_pfc35_analysis.R")
source("../main/step4_subjectselection.R")
source("../main/step5_our_analysis.R")
source("../main/step5_our_analysis_robustness.R")
source("../main/step7_goodness.R")
source("../main/step8_results.R")
source("../main/step9_figures.R")