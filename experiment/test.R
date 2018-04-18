iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

tada <- longitudinalGM::tada
vec <- tada$Gene
vec <- longitudinalGM::symbol_synonyms(vec)
idx <- which(is.na(vec))
tada <- tada[-idx,]; vec <- vec[-idx]
tada$Gene <- vec

#remove duplicated tada by keeping the one with the lowest p-value
tada <- tada[-which(duplicated(tada$Gene)),]

idx <- which(tada$Gene %in% iossifov)
lof <- apply(tada[,c(3:5)], 1, sum)
plot(sort(lof[idx]))

#############
# tool
cutoff_values <- function(report, val = 10:25){
  iossifov <- longitudinalGM::iossifov
  iossifov <- longitudinalGM::symbol_synonyms(iossifov)
  iossifov <- iossifov[!is.na(iossifov)]

  idx <- which(report$Gene %in% iossifov)
  vec <- report$FDR[idx]
  vec <- sort(vec, decreasing = F)
  res <- sapply(val, function(x){
    tmp <- vec[x]
    length(which(report$FDR <= tmp))
  })

  mat <- rbind(val, res, val/res)
  mat[3,] <- sapply(mat[3,], function(x){round(x,3)})
  mat
}
