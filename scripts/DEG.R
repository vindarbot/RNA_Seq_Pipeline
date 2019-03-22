#BiocManager::install("docopt")
library(docopt)

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

library(DESeq2)
library(GOstats)
library(Category)
library(tidyverse)
library(biomaRt)
library(ggplot2)

# On lit le fichier design expérimental, qui a était généré avec Snakemake
xpdesign = read.csv(xpdesign, row.names=1, sep=",")

# On enlève l'extension du nom des fichiers
rownames(xpdesign) <- gsub(".sorted.bam","",rownames(xpdesign))

RPKM=read.csv("RPKM.txt",row.names=1,sep="\t")

colnames(RPKM) <- gsub("Mapping.", "", colnames(RPKM))
colnames(RPKM) <- gsub(".sorted.bam","",colnames(RPKM))

RPKM <- RPKM[,1:ncol(RPKM)-1]

# On lit la matrice de comptage des reads par gène généré par featureCounts
featurescounts=read.csv("featureCounts/counts.txt", sep="", head=T, skip=1, row.names = "Geneid")


#Formatage de la matrice de comptage généré par featureCounts, pour que celle-ci
# puisse être lu par DESeq2 

#Ici on retire les 5 premières colonnes inutiles pour ne garder que les colonnes associées aux nombre 
# de reads par compté pour chaque échantillon
featureMatrix <- featurescounts[ ,6:ncol(featurescounts)]

# Formatage du nom des échantillon 
colnames(featureMatrix) <- gsub("Mapping.", "", colnames(featureMatrix))
colnames(featureMatrix) <- gsub(".sorted.bam","",colnames(featureMatrix))


featureMatrix <- as.matrix(featureMatrix)

# On récupère les conditions expérimentales 
conditionTest <- factor(xpdesign$condition)

# Pour chaque gène, on récupère la valeur médiane du nombre de reads entre les 3 réplicats
# Par exemple, à partir de la matrice de featureCounts, pour le gène AT1G01010:

#           Col_1 Col_2 Col_3 hon4_1 hon4_2 hon4_3
# AT1G01010   201   104   944     95    138    173

# La médiane calculé sera :

#             Col   hon4
#AT1G01010    201    138
MedianeParCondition = t(apply(featureMatrix, 1, tapply, 
                              conditionTest, median)) 

# Ici, on récupère la valeur médiane maximal entre les deux conditions pour chaque gène
# Exemple pour AT1GO1O1O1: 201 (car 201>138)

maxMedian=apply(MedianeParCondition, 1, max) 

# Ainsi, pour chaque gène, si cette valeur médiane maximal est inférieur à 10, on supprime le gène de la matrice de comptage
# Ceci permet d'éliminer les gènes avec très peu de reads, ces gènes 
# pouvant être à tord détectés comme différentiellement exprimés
#
featureMatrix = featureMatrix[maxMedian >= 10,]



# On créé un tableau qui associe pour chaque échantillon le nom de la condition associé (colonne conditionTest)
# Ceci permet de specifier le design expérimental à DESeq2
# Exemple:
#        conditionTest
#Col_2            Col
#Col_1            Col
#Col_3            Col
#hon4_3          hon4
#hon4_1          hon4
#hon4_2          hon4

(coldata <- data.frame(row.names=colnames(featureMatrix[,rownames(xpdesign)]), conditionTest))
as.factor(xpdesign$condition)
colnames(coldata)

# dds permet de transformer la matrice que nous venons de formater en un object DESeq
# ATTENTION à ce que les colonnes de featureMAtrix soit dans le même ordre que le nom
# des lignes (et donc des échantillons) du fichier de design expérimental
# pour s'assurer de ceci, on formate la matrice comme ceci: featureMatrix[,rownames(xpdesign)]

dds <- DESeqDataSetFromMatrix(countData=featureMatrix[,rownames(xpdesign)], colData=coldata, design= ~ conditionTest)

# Cette fonction permet de normaliser le nombre de reads compté selon la taille de la librairie
# de séuqneçage pour chaque échantillon, selon le nombre de total de reads compté
# Ceci permet de normaliser le nombre de reads compté par gène entre les échantillons

dds <- estimateSizeFactors(dds)

dds = estimateDispersions( dds)
dds <- DESeq(dds)
?estimateDispersions
# On transorme le nombre de reads comptés par le log de ce nombre, afin de minimiser les grands écarts
# pouvant être observés entre les gènes avec peu de reads comptés, et les gènes avec beaucoup de reads comptés
rld <- rlog(dds, blind = FALSE)

# On génère une ACP sur les 1000 gènes les + exprimés 
plotPCA(rld, intgroup="conditionTest",ntop = 1000)

# On récupère les résultats obtenus par DESeq2
res <- results(dds)

plotMA(res, ylim=c(-5,5))



## Gènes induits lorsque le gène testé est muté
# Ici on récupère les gènes différentiellement surexprimés dans la condition test,
# ceci selon un seuil de p-value ajustée définie (exemple padj<0.1) et selon un log2FoldChange définit
# (exemple log2FoldCHange de 1, soit un FoldChange de 2)
genes_up = res[ which(res$padj < padj & res$log2FoldChange > lfc), ]

# On fait de même pour les gènes sous-exprimés
genes_down = res[ which(res$padj < padj & res$log2FoldChange < -lfc), ]

genes_up <- genes_up[order(genes_up$padj, decreasing = F),]

## Gènes réprimés lorsque le gène testé est muté
genes_down = res[ which(res$padj < padj & res$log2FoldChange < -lfc), ]

genes_down <- genes_down[order(genes_down$padj, decreasing = F),]

genes_down <- as.data.frame(genes_down)
genes_up <- as.data.frame(genes_up)

# On va récupérer les identifiers TAIR pour chaque gène et les ajouter dans une colonne "ID"
# Ceci sera nécessaire plus tard pour récuperer les Gene Symbols 
identifiants_genes_up <- rownames(genes_up)
identifiants_genes_down <- rownames(genes_down)

genes_up$ID <- identifiants_genes_up
genes_down$ID <- identifiants_genes_down

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)

## Le fichier gene_description a été récupéré sur TAIR10 et nous permet d'accéder aux annotations
# des gènes, à partir de leur identifiants.

# Ici on récupère le fichier de description des gènes , récupéré sur le site de TAIR
gene_description <- read.delim(description,
                               header=FALSE, quote="")

# On formate le nom des identifiants au sein du fichier
genenames = gsub("[.][1234567890]", "",
                 gene_description[,1])

gene_description[,1]=genenames

## Ici on récupère la description des gènes différentiellement exprimés, en regardant les identifants
# des gènes qui "matchent" entre les deux fichiers
genes_match_rows <- match(rownames(genes_signif), gene_description[,1])



# On fait ceci pour les gènes up
genes_match_up <- match(rownames(genes_up), gene_description[,1])

Code_to_description <- gene_description[genes_match_rows,c(1,3)]

Code_to_description_up <- gene_description[genes_match_up,c(1,3)]

colnames(Code_to_description) <- c("Code","Description")

colnames(Code_to_description_up) <- c("Code","Description")

# Et les gènes downs
genes_match_down <- match(rownames(genes_down), gene_description[,1])

Code_to_description_down <- gene_description[genes_match_down,c(1,3)]

colnames(Code_to_description_down) <- c("Code","Description")


# On ajoute la description des gènes dans la variable gene_up et gene_dpwn
genes_up$Description <- Code_to_description_up$Description

genes_down$Description <- Code_to_description_down$Description


# Ici on utilise la librairie bioMart qui permet de traduire et donc de récupérer les Gene Symbols,
# à partir des identifiants TAIR


ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")


TAIRID <- genes_down$ID

# Ici, on déifinit un tableau symbol qui associe chaque TAIR ID avec le GeneSymbol associé, exemple:

#       tair_symbol tair_locus
#1        AT-AER  AT5G16970
#2                AT4G32100
#3                AT2G43120
#4                AT1G30814
#5         PUB29  AT3G18710
#6         APUM6  AT4G25880
symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

# On ajoute dans les variables genes_up et genes_up les genesSymbols associés
# Puis on garde uniquement les colonnes d'interet 
# Enfin, on trie le tableau par padj décroissante:
# Exemple

genes_up <- (merge(x=symbol,y=genes_up,by.x="tair_locus",by.y="ID",all.y=TRUE))
genes_up <- genes_up[c(1,2,4,8,9)]
genes_up <- genes_up[order(genes_up$padj,decreasing = F),]

#     tair_locus tair_symbol log2FoldChange     padj
#25   AT1G36180        ACC2       2.438645 1.252854e-39
#47   AT2G02120      PDF2.1       2.082633 8.393737e-19
#170  AT5G25980        TGG2       1.589085 1.758894e-16
#150  AT4G29190                   1.349187 5.826872e-16

genes_down <- (merge(x=symbol,y=genes_down,by.x="tair_locus",by.y="ID",all.y=TRUE))
genes_down <- genes_down[c(1,2,4,8,9)]
genes_down <- genes_down[order(genes_down$padj,decreasing = F),]

# Enfin on génère le fichiers de sortie des gènes up et down

write_tsv(as.data.frame(genes_up),"DEG/genes_up.txt")
write_tsv(as.data.frame(genes_down),"DEG/genes_down.txt")







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













