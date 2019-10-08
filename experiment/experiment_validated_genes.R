rm(list=ls())
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)

tada <- read.csv("../../raw_data/TADA_Results_2231trios_1333trans_1601cases_5397controls_March26_pvalues.csv") # 18735 genes
tada <- tada[,which(colnames(tada) %in% c("Gene", "dn.LoF", "qvalue", "pval.TADA"))]

validated_genes <- read.csv("../../raw_data/102_genes_20190123.txt", header = F)
validated_genes <- sort(as.vector(validated_genes[,1]))

##################################

length(intersect(colnames(genexp), validated_genes))
length(intersect(tada$Gene, validated_genes))

idx <- which(tada$Gene %in% validated_genes)
