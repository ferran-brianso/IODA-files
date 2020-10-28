# ## ----setup, echo=FALSE, results='hide', message=FALSE, eval=TRUE---------
# oldopt <- options(digits=3)
# options(width=80)
# on.exit( {options(oldopt)} )

## ----instPackages, echo=FALSE, message=FALSE-----------------------------
installifnot <- function (packageName){
  if (!(require(packageName, character.only=TRUE))) {
    install.packages(eval(packageName), dep=TRUE)
  }else{
    print(paste("Package", packageName, "already installed", sep=" "))
  } 
}
installBiocifnot <- function (packageName){
  if (!(require(packageName, character.only=TRUE))) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(eval(packageName), suppressUpdates=TRUE)
  }else{
    print(paste("Package", packageName, "already installed", sep=" "))
  } 
}

## ----directories, echo=FALSE, message=FALSE------------------------------
workingDir <-getwd() # or put here the corresponding working directory
setwd(workingDir)
dataDir<- file.path(workingDir, "dades") 
resultsDir <- file.path(workingDir, "results/mixomics") 

## ----packages, echo=FALSE, message=FALSE, eval=FALSE---------------------
# installifnot("car")
# installifnot("Hmisc")
# installifnot("xtable")
# installifnot("gdata")
# installifnot("gplots")
# installBiocifnot("impute")
# installifnot("mixOmics")
# installifnot("FactoMineR")

## ----loadPreviousData----------------------------------------------------
# Load files
load(file=file.path(dataDir, "raw.data.Rda"), verbose = TRUE)
GenesSelected <- t(gene.data)
ProtsSelected <- t(prot.data)
dim(GenesSelected)
dim(ProtsSelected)

## ----checkPreviousData---------------------------------------------------
# And check that there is none NA value
sum(is.na(GenesSelected))
sum(is.na(ProtsSelected))

# And check if there are some duplicated gene names
geneDups <- rownames(GenesSelected)[which(duplicated(rownames(GenesSelected)))]
if (length(geneDups)>0){
   GenesSelected[which(rownames(GenesSelected) %in% geneDups),]
   GenesSelected <- GenesSelected[!duplicated(rownames(GenesSelected)),]
}
dim(GenesSelected)

# And check if there are some duplicated prot names
protDups <- rownames(ProtsSelected)[which(duplicated(rownames(ProtsSelected)))]
if (length(protDups)>0){
   ProtsSelected[which(rownames(ProtsSelected) %in% protDups),]
   ProtsSelected <- ProtsSelected[!duplicated(rownames(ProtsSelected)),]
}
dim(ProtsSelected)

## ----viewSummary, out.width='.6\\\\linewidth'----------------------------
# View a numeric summary of the data
summary(GenesSelected)
summary(ProtsSelected)
# And a histogram of both data sets
hist(GenesSelected)
hist(ProtsSelected)

## ----scaleData, out.width='.6\\\\linewidth'------------------------------
# Scale all data
genesSC <- scale(GenesSelected)
protsSC <- scale(ProtsSelected) # this is not necessary (prots already centered)
# And view the new summary
summary(genesSC)
summary(protsSC)
# And a histogram of both scaled data sets
hist(genesSC)
hist(protsSC)

## ----transposeData-------------------------------------------------------
# Transpose data to work properly
X <- t(genesSC); dim(X)
Y <- t(protsSC); dim(Y)
head(X)
head(Y)
# And simplify sample names
#rownames(X) <- substr(rownames(X), 3, 9)
#head(X)

## ----plotCorrMatrix------------------------------------------------------
require(mixOmics)
## Show image in R device
#imgCor(X, Y, type='separate', interactive=FALSE)
## Or save it out into a file
png(file=file.path(resultsDir, "corrMatrix.png"), 
    width=800, height=800, res=120)
 imgCor(X, Y)
dev.off()

## ----cvScoreTuning, eval=FALSE, cache=TRUE-------------------------------
grid1 <- seq(0.1, 2, length = 30)
grid2 <- seq(0.001, 0.02, length = 20)
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
## ################################################################
## ## Depending on how fine your grid is, this may take some time
## ################################################################
#cv.score <- tune.rcc(X, Y, grid1=grid1, grid2=grid2)
## ################################################################
###################################################################
###################################################################
#save(cv.score, file=file.path(resultsDir, "cv.score.Rdata"))
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################

## ----loadTunedParams, eval=FALSE-----------------------------------------
## ## This call provided the following values in JAN 2019
## #lambda1 =  0.139
## #lambda2 =  0.0147
## #CVScore =  0.915

## ----rCCA----------------------------------------------------------------
require(mixOmics)
## Run rCCA given the optimal parameters:
if(!(exists("cv.score"))) load(file=file.path(resultsDir, "cv.score.Rdata"))
lambda1 <- cv.score$opt.lambda1
lambda2 <- cv.score$opt.lambda2
result <- rcc(X, Y, ncomp = 3, lambda1 = lambda1, lambda2 = lambda2)
head(result$cor)
head(result$loadings$X)
head(result$loadings$Y)
save(result, file=file.path(resultsDir, "rccResult.Rdata"))


## ----samplesPlot, out.width='.8\\\\linewidth', warning=FALSE-------------
## samples plot
rownames(X)
#duplicated(rownames(X))
#duplicated(colnames(X))
#duplicated(rownames(Y))
#duplicated(colnames(Y))

require(mixOmics)
data("breast.TCGA")
#dim(breast.TCGA$data.train$mrna)

col.groups <- as.numeric(breast.TCGA$data.train$subtype)
col.groups
plotIndiv(result, comp = 1:2, col = col.groups, cex=3)
## Save individual samples plot into a file
png(file=file.path(resultsDir, "indivSamples.png"), 
    width=800, height=800, res=120)
plotIndiv(result, comp = 1:2, col = col.groups, cex=3)
dev.off()

## ----corrCirclePlot05, out.width='.8\\\\linewidth'-----------------------
plotVar(result, comp = 1:2, cutoff=0.5, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.05.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.5, cex = c(3, 3))
dev.off()

## ----corrCirclePlot08, out.width='.8\\\\linewidth'-----------------------
plotVar(result, comp = 1:2, cutoff=0.6, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.06.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.6, cex = c(3, 3))
dev.off()

## ----corrCirclePlot08, out.width='.8\\\\linewidth'-----------------------
plotVar(result, comp = 1:2, cutoff=0.7, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.07.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.7, cex = c(3, 3))
dev.off()

## ----corrCirclePlot08, out.width='.8\\\\linewidth'-----------------------
#plotVar(result, comp = 1:2, cutoff=0.8, cex = c(3, 3))
## Save corr circle plot into a file
#png(file=file.path(resultsDir, "corrCircle.08.png"), 
#    width=800, height=800, res=120)
#plotVar(result, comp = 1:2, cutoff=0.8, cex = c(3, 3))
#dev.off()


## ----relNetworks, out.width='.8\\\\linewidth'----------------------------
#netw.threshold <- 0.8
#network(result, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
### Save relevance network into a file
#png(file=file.path(resultsDir, 
#                   paste("relNetwork", as.character(netw.threshold),
#                         ".png", sep="")),
#    width=800, height=800, res=120)
#network(result, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
#dev.off()


## PREPARE result in order to avoid having duplicated vertex names in the resulting network
## That means adding a  'g' or 'p' to each of the feature names in the rcc result object
colnames(result$X) <- paste0("g.", colnames(result$X))
colnames(result$X)
result$names$colnames$X <- colnames(result$X)

colnames(result$Y) <- paste0("p.", colnames(result$Y))
colnames(result$Y)
result$names$colnames$Y <- colnames(result$Y)



netw.threshold <- 0.5

net <- network(result, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
dev.off()

## Save relevance network into a file
png(file=file.path(resultsDir, 
                   paste("relNetwork", as.character(netw.threshold),
                         ".png", sep="")),
    width=800, height=800, res=120)
network(result, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
dev.off()

## Export it to Cytoscape-readable .graphml format
require(igraph)
write.graph(net$gR, 
            file=file.path(resultsDir, 
                           paste("relNetwork", 
                                 as.character(netw.threshold),".graphml",
                                 sep="")), 
            format = "graphml")

## And obtain and save the specific pairs having corr value over the threshold
ok <- net$M
ok[abs(ok)<netw.threshold] <- NA #drop cases below the threshold
ok <- as.data.frame(as.table(ok)) #turn into a 3-column table
ok <- na.omit(ok)
colnames(ok) <- c("Gene", "Prot", "CorrValue")
ok <- ok[order(-abs(ok$CorrValue)),]    #sort by highest abs corr value
write.csv2(ok, file=file.path(resultsDir, 
                              paste("relNetworkValues", 
                                    as.character(netw.threshold),".csv", 
                                    sep="")),
           row.names = TRUE)




## ----cim-----------------------------------------------------------------
## ----cim, out.width='.8\\linewidth'---------------------------------
cim(result, comp = 1:3, ylab = "genes", xlab = "prots", margins = c(5, 6))
## Save cim heatmap into a file
png(file=file.path(resultsDir, "cimHeatMap.png"), 
    width=800, height=800, res=120)
cim(result, comp = 1:3, ylab = "genes", xlab = "prots", margins = c(5, 6))
dev.off()

