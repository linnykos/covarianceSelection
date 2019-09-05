# goodness of fit
rm(list=ls())
load("../results/step5_res.RData")
diag<- summary_results$diagnostic

red <- rgb(219, 51, 64, max = 255)

pdf("../figures/figure_10.pdf", height = 4.5, width = 4)
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     xlab = "Theoretical quantiles", ylab = "Observed quantiles",
     main = "QQ-plot based on data\nfrom selected partitions")
lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = red)
points(sort(diag), seq(0, 1, length.out = length(diag)), pch = 16)
graphics.off()
