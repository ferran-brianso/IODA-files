library(mixOmics)
data("breast.TCGA")

#View(breast.TCGA)
breast.TCGA$data.train$mirna
# miRNA: data frame with 150 (70) rows and 184 columns in the training (test) data set. 
#        The expression levels of 184 miRNA.

breast.TCGA$data.train$mrna
# mRNA: data frame with 150 (70) rows and 520 columns in the training (test) data set. 
#       The expression levels of 200 mRNA.


breast.TCGA$data.train$protein
# protein: data frame with 150 (70) rows and 142 columns in the training data set only. 
#          The abundance of 142 proteins.


str(breast.TCGA$data.train$mirna)
str(breast.TCGA$data.train$mrna)
str(breast.TCGA$data.train$protein)

genes <- breast.TCGA$data.train$mrna
head(genes)
## samples in rows
## gene ids in columns
length(unique(colnames(genes))) ## 200 genes
colnames(genes)
heatmap(genes)

prots <- breast.TCGA$data.train$protein
head(prots)
## samples in rows (same as in mrna data)
## prot names in columns
length(unique(colnames(prots))) ## 142 genes
colnames(prots)
heatmap(prots)



library(topGO)
BPterms <- ls(GOBPTerm)
head(BPterms)


annot.list <- annFUN.org("BP", feasibleGenes = colnames(genes), 
                         mapping = "org.Hs.eg.db", ID = "symbol")

length(annot.list)

lengths(annot.list) ## the number of genes annotated to each GO term
which(lengths(annot.list)>10)

BPgt10 <- annot.list[which(lengths(annot.list)>10)]
BPgt10 # GO BP terms with >10 genes annotated
names(BPgt10)
head(genes)

## transposo la matriu de gens, pq els gens estiguin en les files
expr.matrix <- t(genes)
dim(expr.matrix)
expr.matrix[1:10,1:8]

num.genes <- length(rownames(expr.matrix))
num.genes

num.samples <- length(colnames(expr.matrix))
num.samples

num.categs <- length(BPgt10)

categ.matrix <- matrix(0, nrow = num.genes, ncol = num.categs)
dim(categ.matrix)
rownames(categ.matrix) <- rownames(expr.matrix)
colnames(categ.matrix) <- names(BPgt10)
categ.matrix[1:20,]

BPgt10

for (j in 1:num.categs){
  categ.matrix[which(rownames(categ.matrix) %in% BPgt10[[j]]), j] <- 1
}

categ.matrix[1:20,]
head(categ.matrix)
colSums(categ.matrix)

annot.matrix <- cbind(expr.matrix, categ.matrix)
head(annot.matrix)
dim(annot.matrix)


###################################################
## function definition: expandAnnotatedMatrix
###################################################
expandAnnotatedMatrix <- function(x, s.cols, c.cols){
  ## s.cols must be higher than dim(x)[2] !!
  for (j in 1:length(c.cols)){
    sums <- apply(x[,s.cols], 2, function(y){sum(y[which(x[,c.cols[j]]==1)])})
    #print(sums)
    if (j==1){
      em <- rbind(x, c(sums,rep(NA, length(c.cols))))
      em[dim(em)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }else{
      em <- rbind(em, c(sums,rep(NA, length(c.cols))))
      em[dim(em)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }
    #rownames(em)[dim(em)[1]] <- sprintf("A%d", j)
    rownames(em)[dim(em)[1]] <- colnames(x)[length(s.cols)+j]
  }
  return(em)
}

###################################################
## function call
###################################################
expand.matrix <- expandAnnotatedMatrix(
  annot.matrix, ## matrix including expr values and binary annotations to bio categs
  1:num.samples, ## range of columns containing sample-gene values (1:N in the demo)
  (num.samples+1):(num.samples+num.categs) ## range of columns containing categs (the rest of the cols)
)

expand.matrix
tail(expand.matrix)
ann.expr.matrix <- expand.matrix[ , 1:num.samples] 
dim(ann.expr.matrix)
tail(ann.expr.matrix)

