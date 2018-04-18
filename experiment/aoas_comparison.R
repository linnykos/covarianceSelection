dat = read.csv("../../../../../Downloads/aoas844_supp.csv", sep = ";")
idx = which(dat$Risk_early == 1)
genes = dat$Gene[idx]

iossifov <- longitudinalGM::iossifov
iossifov <- longitudinalGM::symbol_synonyms(iossifov)
iossifov <- iossifov[!is.na(iossifov)]

length(intersect(iossifov, genes))
