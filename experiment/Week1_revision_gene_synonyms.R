rm(list=ls())

## from https://www.biostars.org/p/14971/
## there seems to be more suggestions in https://support.bioconductor.org/p/74297/
library(org.Hs.eg.db)
# set up your query genes
queryGeneNames <- c('WHRN', 'SANS', 'DFNB31')

# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
result <- aliasSymbol[which(aliasSymbol[,2] %in% queryGeneNames),5]
result

#################3

# let's give a run
load("../../raw_data/newGenexp.RData")
rownames(genexp) <- genexp[,1]
genexp <- genexp[,-1]
genexp <- t(genexp)
genexp <- as.data.frame(genexp)

res <- symbol_synonyms(colnames(genexp), verbose = T)
table(sapply(res, length)) # 347 0's
idx <- which(sapply(res, length) == 0)
res <- res[-idx]
genexp <- genexp[-idx]

# now to deal with non-unique rows
genexp2 <- average_same_columns(genexp)


##############3

# table(sapply(res, length)) #hm, why are some more than 1?
# idx <- which(sapply(res, length) == 6)[1] #  
# res[idx] # "MIR18A"  "MIR19A"  "MIR19B1" "MIR20A"  "MIR92A1" "MIR17HG"
# colnames(genexp)[idx] #"SCAP"
# 
# dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
# sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
# 
# rowidx <- which(aliasSymbol[,2] %in%  colnames(genexp)[idx])
# rowidx2 <- which(aliasSymbol[,5] %in%  colnames(genexp)[idx])
# aliasSymbol[rowidx,]
# aliasSymbol[rowidx2,]
# 
# rowidx <- which(aliasSymbol[,2] %in%  c("FAM46C", "MGEA5", "MKL1", "MKL2", "LHFP", "KIAA1524"))
# aliasSymbol[rowidx,]
# 
# aliasSymbol[which(aliasSymbol[,2] %in% "MIR20A"),]
# 
