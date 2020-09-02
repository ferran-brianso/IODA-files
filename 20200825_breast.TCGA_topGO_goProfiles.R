## carreguem el paquet miXomics i el data set de dades

library(mixOmics)
data("breast.TCGA")
str(breast.TCGA)

## Un cop carregades i revisades les dades, tenim 200 gens i 142 prots 
## per a les mateixes 150 mostres. 

gene.data <- breast.TCGA$data.train$mrna
class(gene.data)
dim(gene.data)

prot.data <- breast.TCGA$data.train$protein
class(prot.data)
dim(prot.data)

## Les proteïnes les hem de linkar primer a gens, per procedir amb l'anotació a GO.
## Això ho faig visitant www.uniprot.org i fent les cerques per humà allà, 
## i he obtingut els corresponents gene IDs:

prot2gene <- read.csv("prots_mapped.txt", header = F, sep = "\t")
prot2gene

## Filtro prot.data perque inclogui sols les prots per a les que tinc gene symbol
prot.data[1:10,1:5]
prot.idx.mapped <- which(colnames(prot.data) %in% prot2gene$V1)
prot.data <- prot.data[ ,prot.idx.mapped]
dim(prot.data)
prot.data[1:10,1:5]

## I els faig canvi de nom perque mostrin els symbols corresponents
colnames(prot.data) <- prot2gene$V2
prot.data[1:10,1:5]
gene.data[1:10,1:5]


## Per afegir info biològica als gene symbols no hi ha problema, 
## ja que podem utilitzar la funció annFUN.org del paquet topGO,
## per anotar els symbols a GO, i quedar-nos llavors aquelles GO:BP que 
## estiguin anotades a 10,15... o més dels gens/prots de les nostres llistes.

library(topGO)

gene.ann <- annFUN.org("BP", feasibleGenes = colnames(gene.data), mapping = "org.Hs.eg.db", ID = "symbol")
gene.ann.over10 <- gene.ann[which(lengths(gene.ann)>=10)]
gene.ann.over10
length(gene.ann.over10)

prot.ann <- annFUN.org("BP", feasibleGenes = colnames(prot.data), mapping = "org.Hs.eg.db", ID = "symbol")
prot.ann.over10 <- prot.ann[which(lengths(prot.ann)>=10)]
prot.ann.over10
length(prot.ann.over10)

## FINALMENT, APLICO EL CRITERI DE QUEDAR-ME CATEGORIES ON LA SUMA D'ANNOTACIONS >= A 15 (sumant gens i prots)
join.ann <- annFUN.org("BP", feasibleGenes = c(colnames(gene.data), colnames(prot.data)), mapping = "org.Hs.eg.db", ID = "symbol")
join.ann.over15 <- join.ann[which(lengths(join.ann)>=15)]
#join.ann.over10 <- c(gene.ann.over10, prot.ann.over10)
#length(join.ann.over10)
#length(unique(names(join.ann.over10)))

## Procedim ara a crear les matrius amb les annotacions binàries

## Primer per les dades d'expressió gènica
gene.matrix <- t(gene.data)

## transposo la matriu de dades, pq els gene symbols estiguin en les files
dim(gene.matrix)
gene.matrix[1:10,1:8]

## determino el nombre de gens, mostres i categories
num.genes <- length(rownames(gene.matrix))
num.genes
num.samples <- length(colnames(gene.matrix))
num.samples
num.gene.categs <- length(join.ann.over15) ## fem servir aqui les categs unides de gens i prots!

## inicialitzo la matriu d'anotacions
gene.categ.matrix <- matrix(0, nrow = num.genes, ncol = num.gene.categs)
dim(gene.categ.matrix)
rownames(gene.categ.matrix) <- rownames(gene.matrix)
colnames(gene.categ.matrix) <- names(join.ann.over15)  ## fem servir aqui les categs unides de gens i prots!
gene.categ.matrix[1:15,]

## omplo les dades de la matriu d'anotacions
for (j in 1:num.gene.categs){
  gene.categ.matrix[which(rownames(gene.categ.matrix) %in% join.ann.over15[[j]]), j] <- 1
}

## verifico com s'ha quedat la matriu d'anotacions
gene.categ.matrix[1:15, ]
colSums(gene.categ.matrix)

## ajunto la matriu d'anotacions a la d'expressió
gene.ann.matrix <- cbind(gene.matrix, gene.categ.matrix)
head(gene.ann.matrix)
dim(gene.ann.matrix)


## Repeteixo el procés per a les dades de proteïnes
prot.matrix <- t(prot.data)

## transposo la matriu de dades, pq els gene symbols estiguin en les files
dim(prot.matrix)
prot.matrix[1:10,1:8]

## determino el nombre de gens, mostres i categories
num.prots <- length(rownames(prot.matrix))
num.prots
num.samples <- length(colnames(prot.matrix))
num.samples
num.prot.categs <- length(join.ann.over15)  ## fem servir aqui les categs unides de gens i prots!

## inicialitzo la matriu d'anotacions
prot.categ.matrix <- matrix(0, nrow = num.prots, ncol = num.prot.categs)
dim(prot.categ.matrix)
rownames(prot.categ.matrix) <- rownames(prot.matrix)
colnames(prot.categ.matrix) <- names(join.ann.over15) ## fem servir aqui les categs unides de gens i prots!
prot.categ.matrix[1:15, 1:8]

## omplo les dades de la matriu d'anotacions
for (j in 1:num.prot.categs){
  prot.categ.matrix[which(rownames(prot.categ.matrix) %in% join.ann.over15[[j]]), j] <- 1
}

## verifico com s'ha quedat la matriu d'anotacions
prot.categ.matrix[1:15, ]
colSums(prot.categ.matrix)

## ajunto la matriu d'anotacions a la d'expressió
prot.ann.matrix <- cbind(prot.matrix, prot.categ.matrix)
head(prot.ann.matrix)
dim(prot.ann.matrix)

### Ara si que tenim el mateix nombre de columnes als dos data sets, corresponents a les 150 mostres
### + 30 categories de GO:BP per les que tenim anotacions abundants de gens i/o prots.
colSums(gene.categ.matrix)
colSums(prot.categ.matrix)


### Per fer la comparativa d'anotacions directament amb goProfiles, 
### necessitem reassignar els gene symbols a Entrez IDs, eliminant NAs retornats 
gene.Entrez <- as.character(mapIds(org.Hs.eg.db, rownames(gene.matrix), 'ENTREZID', 'SYMBOL'))
gene.Entrez <- gene.Entrez[which(!is.na(gene.Entrez))]
prot.Entrez <- as.character(mapIds(org.Hs.eg.db, rownames(prot.matrix), 'ENTREZID', 'SYMBOL'))
prot.Entrez <- prot.Entrez[which(!is.na(prot.Entrez))]

### Procedim a fer la comparació de les categories anotades via data set de gens 
### amb les anotades via data set de prots, amb el goProfiles
library(goProfiles)
gene.BP <- basicProfile(gene.Entrez, onto ="BP", level = 2, orgPackage="org.Hs.eg.db")
prot.BP <- basicProfile(prot.Entrez, onto ="BP", level = 2, orgPackage="org.Hs.eg.db")

printProfiles(gene.prot.BP, percentage = TRUE)
plotProfiles (gene.prot.BP, aTitle="Comparison Genes vs Prots (both as Entrez IDs")
plotProfiles (gene.prot.BP, percentage=T, aTitle="Genes vs Prots", legend=T) 

