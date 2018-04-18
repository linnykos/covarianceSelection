rm(list=ls())
source("../main/step0_header.R")
source("../main/step1_loading.R")

additional_name <- ""
alpha_seq <- seq(0, 0.35, by = 0.025)

for(alpha in alpha_seq){
  print(alpha)
  print(additional_name)
  vars <- ls()
  vars <- vars[-which(vars %in% c("alpha_seq", "alpha", "additional_name"))]
  rm(list = vars)
  load(paste0(save_filepath, "step1_res.RData"))
 
  set.seed(10)
  source("../main/step2_subjectselection.R")
  source("../main/step3_graph.R")
  source("../main/step4_hmrf.R")
  source("../main/step5_results.R")

  save.image(file = paste0(save_filepath, "step5_res_alpha",alpha, additional_name, ".RData"))
}

warnings()
quit(save = "no")
