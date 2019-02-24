library(Rsamtools)
library(GenomicFeatures)
library(DESeq2)
library(GOstats)
library(Category)
library(org.At.tair.db)
library(tidyverse)
library(GenomicAlignments)
library(biomaRt)
library(ggplot2)
library(DEXSeq)
library(docopt)

# library(stringr)
#library(limma)
dir = setwd("~/Desktop/snakemake/Salmon")

BiocManager::install("docopt")


## Pour récupérer la liste des alignements

bamfiles = BamFileList(dir(path = "../Mapping/", pattern = "sorted.bam$"), yieldSize=100000)

## Pour récupérer le fichier annoté des gènes d'Arabidopsis, afin de compter le nombre de reads par gène.
tair10_genes <- makeTxDbFromGFF("reference_AA.gtf", format="gtf", dataSource="TAIR", organism="Arabidopsis thaliana")

## On regroupe selon les gènes.
(ebg = exonsBy(tair10_genes, by="gene"))

## Design expérimental, utile par la suite pour effectuer l'analyse de gènes différentiellement exprimés entre les deux conditions.
xpdesign = read.csv("experimentalDesign.txt", row.names=1, sep=",")



dir = setwd("~/Desktop/SNAKEMAKE")

## La fonction summarizeOverlaps prend le fichiers bam ainsi que le fichier annoté de gènes 
## afin de compter le nombre de reads correspondant à un gène.  Le mode union signifie que si une lecture 
## correspond à une zone où deux gènes se chevauchent, elle n'est pas comptée.
## singleEnd = F dans notre cas (car paired-ends)
## fragments = T car on veut compter les reads d'une pair ayant aligné même si l'autre read de la pair 
## n'a pas été aligné

se = summarizeOverlaps(features=ebg,
                       reads=bamfiles, mode="Union", singleEnd=F, 
                       ignore.strand=TRUE, fragments = T )


dir = setwd("Mapping/")
dir = setwd("../")

## Création d'un tableau du nombre de reads compté
counts = assay(se) 
countgenenames = gsub("[.][1234567890]", "", row.names(counts)) 
rownames(counts)=countgenenames 

lala <- as.data.frame(counts_norm)

sum(lala$Control_7.sorted.bam)

test <- as.matrix(RPKM)
head(test)
pca <- prcomp(t(test))

rownames(pca$x) <- c("SA1","SA2","SA3","SA4","SA5","SA6","SA7","SA8","SA9","SA10","MU1","MU2","MU3")
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + 
  geom_text() +
  xlab(paste("PC1 -", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 -", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")




## Afin d'éliminer les gènes exprimés à très faible niveau, (qui pourraient être prédits à tord
## comme différentiellement exprimés), on retire au sein des deux conditions les gènes dont le 
## comptage médian est inférieur à 10.

MedianeParCondition = t(apply(featureMatrix, 1, tapply, 
                              xpdesign, median)) 

maxMedian=apply(MedianeParCondition, 1, max) 

featureMatrix = featureMatrix[maxMedian >= 10,]











### FEATURE COUNTS

featurescounts=read.csv("counts.txt", sep="", head=T, skip=1, row.names = "Geneid")



# RPKM=read.csv("~/Desktop/Analyse_Col0/RPKM.txt",row.names=1,sep="\t")

# RPKM <- RPKM[,1:ncol(RPKM)-1]


## PCA 

test <- as.matrix(RPKM)
head(test)
pca <- prcomp(t(test))

rownames(pca$x) <- c("SA1","SA2","SA3","SA4","SA8","SA9","SA10","MU1","MU2","MU3")
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + 
  geom_text() +
  xlab(paste("PC1 -", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 -", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")



### DEG


featureMatrix <- featurescounts[ ,6:ncol(featurescounts)]
  
colnames(featureMatrix) <- gsub("Mapping.", "", colnames(featureMatrix))

featureMatrix <- as.matrix(featureMatrix)

conditionTest <- factor(c(rep("Col",2), rep("Mutant", 3)))

MedianeParCondition = t(apply(test, 1, tapply, 
                              conditionTest, median)) 

maxMedian=apply(MedianeParCondition, 1, max) 

featureMatrix = featureMatrix[maxMedian >= 10,]




(coldata <- data.frame(row.names=colnames(featureMatrix), conditionTest))

dds = DESeqDataSetFromMatrix(countData=featureMatrix, colData=coldata, design= ~ conditionTest)

dds <- DESeq(dds)

rld <- rlog(dds, blind = FALSE)


head(assay(rld))
hist(assay(rld))


plotPCA(rld, intgroup="conditionTest",ntop = 1000)


res <- results(dds)

res <- results(dds, contrast=c("conditionTest","col0_kaku","mutant"))
res <- results(dds, contrast=c("conditionTest","col0_hon","col0_kaku"))
res.ctr2.ctr3 <- results(dds, contrast=c("conditionTest","ctr2","ctr3"))









## Gènes différentiellement exprimés

## On se sert du package DESeq2 afin d'analyser les gènes différentiellement (test de student entre les deux conditions)
# Dans notre cas, 2 conditions (contrôle et mutant)
# Afin d'avoir plus de robustesse au niveau du test statistique, j'ai décidé de regrouper SA1 SA2 SA3 SA4 comme contrôles
# (et non traiter différement SA1 et SA2 d'une part et SA3 et SA4 d'autre part)

# Création d'une matrice utilisable par DESeq2 pour la suite. On spécifie les données ainsi que 
# le design expérimental
cds = DESeqDataSetFromMatrix(countData=counts,
                             colData=xpdesign,
                             design= ~ condition)



## NORMALISATION: La fonction estimateSizeFactors permet de normaliser le nombre de reads en fonction 
# du nombre total de reads par échantillon

cds = estimateSizeFactors(dds)

cds = estimateDispersions( dds )
colData(cds)
plotDispEsts( cds )
?DESeqDataSetFromMatrix
## On lance l'analyse
cds <- DESeq(cds)

## On garde les résultats dans une variable
res <- results(dds)


res <- results(dds, contrast=c("conditionTest","ctr2","exp"))
res <- results(dds, contrast=c("conditionTest","ctl","ctr3"))
res <- results(dds, contrast=c("conditionTest","ctr2","ctr3"))

#########
rld <- rlog( cds )
rld <- rlogTransformation(cds, blind=TRUE)



as.data.frame(rld)

sampleDists <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                      rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
sum(res$padj < 0.05, na.rm=T)



plotPCA(rld,ntop=500,intgroup="condition")

rld

plotMA(res, ylim=c(-5,5))
?plotPCA



## Gènes induits lorsque le gène testé est muté
genes_up = res[ which(res$padj < 0.11 & res$log2FoldChange > 1), ]

genes_up <- genes_up[order(genes_up$padj, decreasing = F),]

## Gènes réprimés lorsque le gène testé est muté
genes_down = res[ which(res$padj < 0.11 & res$log2FoldChange < -1), ]

genes_down <- genes_down[order(genes_down$padj, decreasing = F),]

genes_down <- as.data.frame(genes_down)
genes_up <- as.data.frame(genes_up)

identifiants_genes_up <- rownames(genes_up)
identifiants_genes_down <- rownames(genes_down)
identifiants <- c(identifiants_genes_down,identifiants_genes_up)

genes_up$ID <- identifiants_genes_up
genes_down$ID <- identifiants_genes_down

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)


## Le fichier gene_description a été récupéré sur TAIR10 et nous permet d'accéder aux annotations
# des gènes, à partir de leur identifiants.

gene_description <- read.delim("~/Desktop/RNA_SEQ/gene_description_20131231.txt",
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

genes_up <- (merge(x=symbol,y=genes_up,by.x="tair_locus",by.y="ID"))
genes_up <- genes_up[c(1,2,4,8,9)]
genes_up <- genes_up[order(genes_up$padj,decreasing = F),]

genes_down <- (merge(x=symbol,y=genes_down,by.x="tair_locus",by.y="ID"))
genes_down <- genes_down[c(1,2,4,8,9)]
genes_down <- genes_down[order(genes_down$padj,decreasing = F),]

  
write_tsv(as.data.frame(genes_up),"~/Desktop/genes_up.txt")
write_tsv(as.data.frame(genes_down),"~/Desktop/genes_down.txt")
write_tsv(as.data.frame(identifiants),"~/Desktop/id.txt")







### Analyse d'enrichissement

params=new("GOHyperGParams",
           geneIds=rownames(genes_signif),
           universeGeneIds=rownames(counts_filtered),
           annotation="org.At.tair",
           ontology="MF",
           pvalueCutoff=0.01,
           conditional=TRUE,
           testDirection="over")

(overRepresented=hyperGTest(params))

summary(overRepresented)[,c(1,2,5,6,7)]



counts_norm <- counts(cds, normalized = T)












?DEXSeq

