B = matrix(c(.25, .5, .25, .5, .25, .25, .25, .25, 1/6), 3, 3)
zz <- eigen(B)
zz$vector[,1] * sqrt(zz$value[1])
zz$vector[,3] * sqrt(abs(zz$value[3]))

.61^2-.35^2 # so to evaluate probability, you use POS - NEG, and to evaluate distance, you just use Euclidean distance


