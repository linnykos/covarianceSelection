rm(list=ls())
devtools::install_github("linnylin92/longitudinalGM",
                         subdir = "longitudinalGM", force = T)

source("../experiment/simulation_base.R")
source("../experiment/simulation_helper.R")

library(longitudinalGM)
set.seed(10)
