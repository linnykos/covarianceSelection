rm(list=ls())
alpha <- 0.1

additional_name_vec <- c("", "_alternative", "_pfc35")

for(simult_idx in 1:length(additional_name_vec)){
  print(paste0("===Master iteration ", simult_idx, "==="))
  var <- ls()
  var <- var[-which(var %in% c("simult_idx", "additional_name_vec", "alpha"))]
  rm(list = var)

  additional_name <- additional_name_vec[simult_idx]

  source("../main/step0_header.R")
  warnings()
  source("../main/step1_loading.R")
  warnings()

  set.seed(10)
  source("../main/step2_subjectselection.R")
  warnings()
  source("../main/step3_graph.R")
  warnings()
  source("../main/step4_hmrf.R")
  warnings()
  source("../main/step5_results.R")
  warnings()
}

warnings()
quit(save = "no")
