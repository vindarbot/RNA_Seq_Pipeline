library(DESeq2)
library(GOstats)
library(Category)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(docopt)

#BiocManager::install("docopt")

"Analyse des gènes différentiellement exprimés (DESeq2)
Usage:
  DEG.R --lfc=<lfc> --padj=<padj> --xpdesign=<xpdesign> --description=<description> --outprefix=<outprefix> 
  Options:
    --lfc=<lfc>                      Log 2 Fold change à appliquer entre les deux contions
    --padj=<padj>                    Seuil de pvalue ajustée à appliquer 
    --xpdesign=<xpdesign>            Design expérimental de l'expérience
    --description=<description>      Fichier d'annoation des gènes
    --outprefix=<outprefix>          Dossier avec les résultats
" -> doc

opts <- docopt(doc)

lfc <- opts[['lfc']]
lfc <- as.numeric(lfc)  
padj <- opts[['padj']]
xpdesign <- opts[['xpdesign']]
outprefix <- opts[['outprefix']]
description <- opts[['description']]


xpdesign = read.csv(xpdesign, row.names=1, sep=",")

featurescounts=read.csv("featureCounts/counts.txt", sep="", head=T, skip=1, row.names = "Geneid")


### Formatage de la matrice de comptage généré par featureCounts, pour que celle-ci
# puisse être lu par DESeq2 

featureMatrix <- featurescounts[ ,6:ncol(featurescounts)]
  
colnames(featureMatrix) <- gsub("Mapping.", "", colnames(featureMatrix))

featureMatrix <- as.matrix(featureMatrix)

conditionTest <- factor(xpdesign$condition)

MedianeParCondition = t(apply(featureMatrix, 1, tapply, 
                              conditionTest, median)) 

maxMedian=apply(MedianeParCondition, 1, max) 

featureMatrix = featureMatrix[maxMedian >= 10,]




(coldata <- data.frame(row.names=colnames(featureMatrix[,rownames(xpdesign)]), conditionTest))
as.factor(xpdesign$condition)
colnames(coldata)

dds <- DESeqDataSetFromMatrix(countData=featureMatrix[,rownames(xpdesign)], colData=coldata, design= ~ conditionTest)

dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)

dds = estimateDispersions( dds )

rld <- rlog(dds, blind = FALSE)


plotPCA(rld, intgroup="conditionTest",ntop = 1000)

res <- results(dds)

plotMA(res, ylim=c(-5,5))



## Gènes induits lorsque le gène testé est muté
genes_up = res[ which(res$padj < padj & res$log2FoldChange > lfc), ]

genes_up <- genes_up[order(genes_up$padj, decreasing = F),]

## Gènes réprimés lorsque le gène testé est muté
genes_down = res[ which(res$padj < padj & res$log2FoldChange < -lfc), ]

genes_down <- genes_down[order(genes_down$padj, decreasing = F),]

genes_down <- as.data.frame(genes_down)
genes_up <- as.data.frame(genes_up)

identifiants_genes_up <- rownames(genes_up)
identifiants_genes_down <- rownames(genes_down)

genes_up$ID <- identifiants_genes_up
genes_down$ID <- identifiants_genes_down

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)

## Le fichier gene_description a été récupéré sur TAIR10 et nous permet d'accéder aux annotations
# des gènes, à partir de leur identifiants.

gene_description <- read.delim(description,
                               header=FALSE, quote="")

genenames = gsub("[.][1234567890]", "",
                 gene_description[,1])

gene_description[,1]=genenames

## Ici on récupère la description des gènes différentiellement exprimés
genes_match_rows <- match(rownames(genes_signif), gene_description[,1])

genes_match_up <- match(rownames(genes_up), gene_description[,1])

Code_to_description <- gene_description[genes_match_rows,c(1,3)]

Code_to_description_up <- gene_description[genes_match_up,c(1,3)]


colnames(Code_to_description) <- c("Code","Description")

colnames(Code_to_description_up) <- c("Code","Description")


genes_match_down <- match(rownames(genes_down), gene_description[,1])

Code_to_description_down <- gene_description[genes_match_down,c(1,3)]

colnames(Code_to_description_down) <- c("Code","Description")


genes_up$Description <- Code_to_description_up$Description

genes_down$Description <- Code_to_description_down$Description


# Pour récupérer les ID symbols à partir des TAIR ID.

ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")

TAIRID <- genes_down$ID

symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

genes_up <- (merge(x=symbol,y=genes_up,by.x="tair_locus",by.y="ID",all.y=TRUE))
genes_up <- genes_up[c(1,2,4,8,9)]
genes_up <- genes_up[order(genes_up$padj,decreasing = F),]

genes_down <- (merge(x=symbol,y=genes_down,by.x="tair_locus",by.y="ID",all.y=TRUE))
genes_down <- genes_down[c(1,2,4,8,9)]
genes_down <- genes_down[order(genes_down$padj,decreasing = F),]


write_tsv(as.data.frame(genes_up),outprefix+"genes_up.txt")
write_tsv(as.data.frame(genes_down),outprefix+"genes_down.txt")







# ### Analyse d'enrichissement

# params=new("GOHyperGParams",
#            geneIds=rownames(genes_signif),
#            universeGeneIds=rownames(counts_filtered),
#            annotation="org.At.tair",
#            ontology="MF",
#            pvalueCutoff=0.01,
#            conditional=TRUE,
#            testDirection="over")

# (overRepresented=hyperGTest(params))

# summary(overRepresented)[,c(1,2,5,6,7)]



# counts_norm <- counts(cds, normalized = T)













