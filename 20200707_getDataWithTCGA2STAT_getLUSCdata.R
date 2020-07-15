## Get data from TCGA
lusc.rnaseq <- getTCGA(disease="LUSC", data.type="RNASeq", type="RPKM", clinical=TRUE)
lusc.mirna <- getTCGA(disease="LUSC", data.type="miRNASeq", type="rpmmm", clinical=TRUE)

## Select samples with matched tumor and normal data
lusc.rnaseq.tum.norm <- TumorNormalMatch(lusc.rnaseq$dat)
lusc.mirna.tum.norm <- TumorNormalMatch(lusc.mirna$dat)

colnames(lusc.rnaseq.tum.norm$primary.tumor)
colnames(lusc.mirna.tum.norm$primary.tumor)
which(colnames(lusc.rnaseq.tum.norm$primary.tumor) %in% colnames(lusc.mirna.tum.norm$primary.tumor))

colnames(lusc.rnaseq.tum.norm$normal)
colnames(lusc.mirna.tum.norm$normal)
which(colnames(lusc.rnaseq.tum.norm$normal) %in% colnames(lusc.mirna.tum.norm$normal))

which.rnaseq.tum <- which(colnames(lusc.rnaseq.tum.norm$primary.tumor) %in% colnames(lusc.mirna.tum.norm$primary.tumor))
which.rnaseq.tum
which.rnaseq.norm <- which(colnames(lusc.rnaseq.tum.norm$normal) %in% colnames(lusc.mirna.tum.norm$normal))
which.rnaseq.norm
which.mirna.tum <- which(colnames(lusc.mirna.tum.norm$primary.tumor) %in% colnames(lusc.rnaseq.tum.norm$primary.tumor))
which.mirna.tum
which.mirna.norm <- which(colnames(lusc.mirna.tum.norm$normal) %in% colnames(lusc.rnaseq.tum.norm$normal))
which.mirna.norm

cols.rnaseq.t <- gsub("TCGA", "T", colnames(lusc.rnaseq.tum.norm$primary.tumor)[which.rnaseq.tum])
cols.rnaseq.t
cols.rnaseq.n <- gsub("TCGA", "N", colnames(lusc.rnaseq.tum.norm$normal)[which.rnaseq.norm])
cols.rnaseq.n
#cols.mirna.t <- gsub("TCGA", "T", colnames(lusc.mirna.tum.norm$primary.tumor)[which.mirna.tum])
#cols.mirna.t ## should be the same!

lusc.rnaseq.common <- cbind(lusc.rnaseq.tum.norm$primary.tumor[ ,which.rnaseq.tum], 
                     lusc.rnaseq.tum.norm$normal[ ,which.rnaseq.norm])
colnames(lusc.rnaseq.common)
colnames(lusc.rnaseq.common) <- c(cols.rnaseq.t, cols.rnaseq.n)
colnames(lusc.rnaseq.common)

lusc.mirna.common <- cbind(lusc.mirna.tum.norm$primary.tumor[,which.mirna.tum], 
                            lusc.mirna.tum.norm$normal[,which.mirna.norm])
colnames(lusc.mirna.common)
colnames(lusc.mirna.common) <- colnames(lusc.rnaseq.common) ## should be the same!
colnames(lusc.mirna.common)

class(lusc.rnaseq.common)
dim(lusc.rnaseq.common)
dim(lusc.mirna.common)
head(lusc.rnaseq.common)
head(lusc.mirna.common)

## To reduce computation time, we only include features that are at the top 1% of variance across patients:
## By now, not considering the paired condition of the samples
rnaseq.var <- apply(lusc.rnaseq.common, 1, var)
rnaseq.data <- subset(lusc.rnaseq.common, rnaseq.var >= quantile(rnaseq.var, 0.99, na.rm=T) & !is.na(rnaseq.var))
dim(rnaseq.data)

mirna.var <- apply(lusc.mirna.common, 1, var)
mirna.data <- subset(lusc.mirna.common, mirna.var >= quantile(mirna.var, 0.90, na.rm=T) & !is.na(mirna.var))
dim(mirna.data)

save(file = "lusc.data.topvar.Rda", rnaseq.data, mirna.data)
#load("lusc.data.topvar.Rda", verbose = T)



