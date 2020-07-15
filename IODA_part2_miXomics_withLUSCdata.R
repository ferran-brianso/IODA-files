
## ----directories, echo=FALSE, message=FALSE------------------------------
workingDir <-getwd() # or put here the corresponding working directory
setwd(workingDir)
dataDir<- file.path(workingDir, "dades") 
resultsDir <- file.path(workingDir, "results/mixomics") 

## ----loadPreviousData----------------------------------------------------
# Load files
load(file="lusc.data.topvar.Rda", verbose = TRUE)
dim(rnaseq.data)
dim(mirna.data)
omics1 <- rnaseq.data
omics2 <- mirna.data

## ----checkPreviousData---------------------------------------------------
# And check that there is none NA value
sum(is.na(omics1))
sum(is.na(omics2))

# And check if there are some duplicated row names in 1st matrix
rowDups1 <- rownames(omics1)[which(duplicated(rownames(omics1)))]
if (length(rowDups1)>0){
   omics1[which(rownames(omics1) %in% rowDups1),]
   omics1 <- omics1[!duplicated(rownames(omics1)),]
}
dim(omics1)

# And check if there are some duplicated row names in 2nd matrix
rowDups2 <- rownames(omics2)[which(duplicated(rownames(omics2)))]
if (length(rowDups2)>0){
   omics2[which(rownames(omics2) %in% rowDups2),]
   omics2 <- omics2[!duplicated(rownames(omics2)),]
}
dim(omics2)

## ----viewSummary, out.width='.6\\\\linewidth'----------------------------
# View a numeric summary of the data
summary(omics1)
summary(omics2)
# And a histogram of both data sets
hist(omics1)
hist(omics2)

## ----scaleData, out.width='.6\\\\linewidth'------------------------------
# Scale all data
omics1SC <- scale(omics1)
omics2SC <- scale(omics2) # this is not necessary (prots already centered)
# And view the new summary
summary(omics1SC)
summary(omics2SC)
# And a histogram of both scaled data sets
hist(omics1SC)
hist(omics2SC)

## ----transposeData-------------------------------------------------------
# Transpose data to work properly
X <- t(omics1SC); dim(X)
Y <- t(omics2SC); dim(Y)
head(X)
head(Y)
# And simplify sample names
#rownames(X) <- substr(rownames(X), 3, 9)
#head(X)

## ----plotCorrMatrix------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics", update = F)
require(mixOmics)
## Show image in R device
#imgCor(X, Y, type='separate', interactive=FALSE)
## Or save it out into a file
png(file=file.path(resultsDir, "corrMatrix.png"), 
    width=800, height=800, res=120)
 imgCor(X, Y)
dev.off()

## ----cvScoreTuning, eval=FALSE, cache=TRUE-------------------------------
grid1 <- seq(0.01, 1, length = 50)
grid2 <- seq(0.01, 0.1, length = 20)
## ################################################################
## ## Depending on how fine your grid is, this may take some time
## ################################################################
cv.score <- tune.rcc(X, Y, grid1=grid1, grid2=grid2)
cv.score
## ################################################################
save(cv.score, file=file.path(resultsDir, "cv.score.Rdata"))

## ----loadTunedParams, eval=FALSE-----------------------------------------
## ## This call provided the following values in JUL 2020
## #lambda1 =  0.033
## #lambda2 =  0.081
## #CVScore =  0.565

## ----rCCA----------------------------------------------------------------
## Run rCCA given the optimal parameters:
if(!(exists("cv.score"))) load(file=file.path(resultsDir, "cv.score.Rdata"))
lambda1 <- cv.score$opt.lambda1
lambda2 <- cv.score$opt.lambda2
result <- rcc(X, Y, ncomp = 3, lambda1 = lambda1, lambda2 = lambda2)
head(result$cor)

## ----samplesPlot, out.width='.8\\\\linewidth', warning=FALSE-------------
## samples plot
rownames(X)
col.groups <- rep(1:2, each=9)
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
plotVar(result, comp = 1:2, cutoff=0.8, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.08.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.8, cex = c(3, 3))
dev.off()

## ----corrCirclePlot08, out.width='.8\\\\linewidth'-----------------------
plotVar(result, comp = 1:2, cutoff=0.85, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.085.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.85, cex = c(3, 3))
dev.off()

## ----corrCirclePlot08, out.width='.8\\\\linewidth'-----------------------
plotVar(result, comp = 1:2, cutoff=0.9, cex = c(3, 3))
## Save corr circle plot into a file
png(file=file.path(resultsDir, "corrCircle.09.png"), 
    width=800, height=800, res=120)
plotVar(result, comp = 1:2, cutoff=0.9, cex = c(3, 3))
dev.off()

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

netw.threshold <- 0.65

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
cim(result, comp = 1:3, ylab = "rnaseq", xlab = "mirna", margins = c(5, 6))
## Save cim heatmap into a file
png(file=file.path(resultsDir, "cimHeatMap.png"), 
    width=800, height=800, res=120)
cim(result, comp = 1:3, ylab = "rnaseq", xlab = "mirna", margins = c(5, 6))
dev.off()


### Sparse Projected Latent Structures (sPLS method)
#source: https://www.bioconductor.org/packages/release/bioc/vignettes/mixOmics/inst/doc/vignette.html#pls

dim(X)
dim(Y)

MyResult.spls <- spls(X, Y, keepX = c(30,30), keepY = c(15,15))  

#plotIndiv(MyResult.spls)                                      
plotIndiv(MyResult.spls, group = c(rep("Tumor",9), rep("Normal", 9)),
          rep.space = "XY-variate", legend = TRUE,
          legend.title = 'Sample type',
          title = 'LUSC: sPLS')

plotVar(MyResult.spls, cex=c(3,2), legend = TRUE)
coordinates <- plotVar(MyResult.spls, plot = FALSE)
coordinates

cim(MyResult.spls, comp = 1)

network(MyResult.spls, comp = 1)


plotArrow(MyResult.spls, group = c(rep("Tumor",9), rep("Normal", 9)), legend = TRUE,
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2')


## Block sparse PLS-DA
group = c(rep("Tumor",9), rep("Normal", 9))
group


X <- list(rnaSeq = t(omics1SC),
          miRNA = t(omics2SC))
dim(X$rnaSeq)
dim(X$miRNA)

Y <- as.factor(group)
summary(Y)

list.keepX <- list(rnaSeq = c(20, 15), miRNA = c(15,5))


MyResult.block.splsda <- block.splsda(X, Y, keepX=list.keepX)

circosPlot(MyResult.block.splsda, cutoff = 0.75)

