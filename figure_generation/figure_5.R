rm(list=ls())
source("../simulation/simulation_helper.R")
load("../results/simulation_SnR.RData")

generate_covariance <- function(d, percentage){
  covar_base <- .generate_block(d)
  covar_shuffle1 <- .shuffle(covar_base, percentage = percentage)
  covar_shuffle2 <- .shuffle(covar_base, percentage = percentage)

  list(covar_base = covar_base, covar_shuffle1 = covar_shuffle1,
       covar_shuffle2 = covar_shuffle2)
}

.clockwise90 = function(a) { t(a[nrow(a):1,]) }

set.seed(10)
covar_list1 <- generate_covariance(d = paramMat[1,"d"],
                                  percentage = 0.25)
covar_list2 <- generate_covariance(d = paramMat[1,"d"],
                                   percentage = 0.75)

cov1 <- covar_list1$covar_base
cov2 <- covar_list1$covar_shuffle1
cov3 <- covar_list2$covar_shuffle1

colfunc <- colorRampPalette(c(rgb(247, 234, 200, max = 255),
                              rgb(205,40,54,maxColorValue=255)))
col_vec <- colfunc(100)

pdf("../figures/figure_5.pdf", height = 2.5, width = 7)
par(mfrow = c(1,3), mar = c(0.1, 0.1, 4, 0.1))
image(.clockwise90(cov1), col = col_vec, asp = T, main = "Base Covariance",
      xlab = "", ylab = "", xaxt = "n", bty = "n", yaxt = "n")

image(.clockwise90(cov2), col = col_vec, asp = T, main = "Covariance Shuffled 25%",
      xlab = "", ylab = "", xaxt = "n", bty = "n", yaxt = "n")

image(.clockwise90(cov3), col = col_vec, asp = T, main = "Covariance Shuffled 75%",
      xlab = "", ylab = "", xaxt = "n", bty = "n", yaxt = "n")
graphics.off()
