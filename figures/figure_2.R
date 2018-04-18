# goodness of fit
rm(list=ls())
load("../results/step5_res_alternative.RData")
diag_alternative <- summary_results$diagnostic

load("../results/step5_res_pfc35.RData")
diag_pfc35 <- summary_results$diagnostic

red <- rgb(219, 51, 64, max = 255)

pdf("../figures/figure_2.pdf", height = 4, width = 7)
par(mfrow = c(1,2))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     main = "QQ-plot based on partitions\nfrom only Window 1B",
     xlab = "Theoretical quantiles", ylab = "Observed quantiles")
lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = red)
points(sort(diag_pfc35), seq(0, 1, length.out = length(diag_pfc35)), pch = 16)

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     main = "QQ-plot based on all partitions",
     xlab = "Theoretical quantiles", ylab = "Observed quantiles")
lines(c(0,1), c(0,1), lwd = 2, lty = 2, col = red)
points(sort(diag_alternative), seq(0, 1, length.out = length(diag_alternative)), pch = 16)
graphics.off()
