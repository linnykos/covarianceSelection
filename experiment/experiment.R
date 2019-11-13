rm(list=ls())
load("../results/tmp.RData")
dat_list = dat[1:15]
split1 = 1:7
dat1 <- do.call(rbind, dat_list[split1])
dat2 <- do.call(rbind, dat_list[-split1])

cov1 <- cov(dat1)
cov2 <- cov(dat2)

.clockwise90 = function(a) { t(a[nrow(a):1,]) }

png("../figures/tmp.png", height = 1300, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(1,1,1,1))
image(.clockwise90(cov1))
image(.clockwise90(cov2))
graphics.off()