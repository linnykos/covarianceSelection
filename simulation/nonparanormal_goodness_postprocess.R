rm(list=ls())
load("../results/nonparanormal_goodness.RData")

plot(sort(goodness_our), seq(0,1,length.out = length(goodness_our)), asp = T, xlab = "Theoretical quantiles",
     pch = 16,
     ylab = "Observed quantiles", main = "QQ-plot based on data\nfrom selected partitions")
lines(c(0,1),c(0,1), col = "red", lty = 2)
