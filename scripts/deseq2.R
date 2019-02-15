library(Rsamtools)
library(GenomicFeatures)
library(DESeq2)
library(GOstats)
library(Category)
library(org.At.tair.db)
library(tidyverse)
library(GenomicAlignments)
library(gplots)

# library(stringr)
#library(limma)
dir = setwd("~/Desktop/RNA_SEQ")

# BiocManager::install("pheatmap")


## Pour récupérer la liste des alignements
bamfiles = BamFileList(dir(path = "Alignement/", pattern = "tair10.sorted.bam$"), yieldSize=100000)

## Pour récupérer le fichier annoté des gènes d'Arabidopsis, afin de compter le nombre de reads par gène.
tair10_genes <- makeTxDbFromGFF("TAIR10.gtf", format="gtf", dataSource="TAIR", organism="Arabidopsis thaliana")

## On regroupe selon les gènes.
(ebg = exonsBy(tair10_genes, by="gene"))

## Design expérimental, utile par la suite pour effectuer l'analyse de gènes différentiellement exprimés entre les deux conditions.
xpdesign = read.csv("experimentalDesign.txt", row.names=1, sep=",")

dir = setwd("Alignement")

## La fonction summarizeOverlaps prend le fichiers bam ainsi que le fichier annoté de gènes 
## afin de compter le nombre de reads correspondant à un gène.  Le mode union signifie que si une lecture 
## correspond à une zone où deux gènes se chevauchent, elle n'est pas comptée.
## singleEnd = F dans notre cas (car paired-ends)
## fragments = T car on veut compter les reads d'une pair ayant aligné même si l'autre read de la pair 
## n'a pas été aligné

se = summarizeOverlaps(features=ebg,
                       reads=bamfiles, mode="Union", singleEnd=F, 
                       ignore.strand=TRUE, fragments = T )


dir = setwd("../")

## Création d'un tableau du nombre de reads compté
counts = assay(se) 
countgenenames = gsub("[.][1234567890]", "", row.names(counts)) 
rownames(counts)=countgenenames 

## Afin d'éliminer les gènes exprimés à très faible niveau, (qui pourraient être prédits à tord
## comme différentiellement exprimés), on retire au sein des deux conditions les gènes dont le 
## comptage médian est inférieur à 10.

MedianeParCondition = t(apply(counts, 1, tapply, 
                             xpdesign, median)) 

maxMedian=apply(MedianeParCondition, 1, max) 

counts_filtered = counts[maxMedian >= 10,]



## Gènes différentiellement exprimés

## On se sert du package DESeq2 afin d'analyser les gènes différentiellement (test de student entre les deux conditions)
# Dans notre cas, 2 conditions (contrôle et mutant)
# Afin d'avoir plus de robustesse au niveau du test statistique, j'ai décidé de regrouper SA1 SA2 SA3 SA4 comme contrôles
# (et non traiter différement SA1 et SA2 d'une part et SA3 et SA4 d'autre part)

# Création d'une matrice utilisable par DESeq2 pour la suite. On spécifie les données ainsi que 
# le design expérimental
cds = DESeqDataSetFromMatrix(countData=counts_filtered,
                             colData=xpdesign,
                             design= ~ condition)

## NORMALISATION: La fonction estimateSizeFactors permet de normaliser le nombre de reads en fonction 
# du nombre total de reads par échantillon

cds = estimateSizeFactors(cds)

cds = estimateDispersions( cds )
colData(cds)
plotDispEsts( cds )

## On lance l'analyse
cds <- DESeq(cds)

## On garde les résultats dans une variable
res <- results(cds)


# res.col.con t <- results(cds, contrast=c("condition","col","control"))
# res.col.traitemen t <- results(cds, contrast=c("condition","col","traitement"))
# res.control.traitemen t <- results(cds, contrast=c("condition","control","traitement"))

#########
rld <- rlog( cds )

sampleDists <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
# sum(res$padj < 0.05, na.rm=T)


?results
plotPCA(cds)
plotMA(res, ylim=c(-5,5))

## Gènes induits lorsque le gène testé est muté
genes_up = res[ which(res$padj < 0.05 & res$log2FoldChange > 0), ]


## Gènes réprimés lorsque le gène testé est muté&
genes_down = res[ which(res$padj < 0.05 & res$log2FoldChange < 0), ]

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)

identifiants_genes_up <- rownames(genes_up)
identifiants_genes_down <- rownames(genes_down)
identifiants <- rownames(genes_signif)


write_tsv(as.data.frame(identifiants),"../identifiants.control.traitement.txt")
### Gene Annotation

## Le fichier gene_description a été récupéré sur TAIR10 et nous permet d'accéder aux annotations
# des gènes, à partir de leur identifiants.

gene_description <- read.delim("gene_description_20131231.txt",
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

write_tsv(as.data.frame(Code_to_description_down),"Data/Statistiques/genes_sous_régulés.txt")
write_tsv(as.data.frame(Code_to_description_up),"Data/Statistiques/genes_sur_régulés.txt")

### Analyse d'enrichissement

params=new("GOHyperGParams",
           geneIds=rownames(genes_signif),
           universeGeneIds=rownames(counts_filtered),
           annotation="org.At.tair",
           ontology="MF",
           pvalueCutoff=0.01,
           conditional=TRUE,
           testDirection="over")

params_up=new("GOHyperGParams",
           geneIds=rownames(genes_up),
           universeGeneIds=rownames(counts_filtered),
           annotation="org.At.tair",
           ontology="MF",
           pvalueCutoff=0.01,
           conditional=TRUE,
           testDirection="over")

params_down=new("GOHyperGParams",
           geneIds=rownames(genes_down),
           universeGeneIds=rownames(counts_filtered),
           annotation="org.At.tair",
           ontology="MF",
           pvalueCutoff=0.01,
           conditional=TRUE,
           testDirection="over")

(overRepresented=hyperGTest(params))

(overRepresented_up=hyperGTest(params_up))

(overRepresented_down=hyperGTest(params_down))

summary(overRepresented)[,c(1,2,5,6,7)]

summary(overRepresented_up)[,c(1,2,5,6,7)]




down <- summary(overRepresented_down)[,c(1,2,5,6,7)]
down <- down[order(-down$Count),]
down$Term <- factor(down$Term,levels=down$Term[order(down$Count)])

graph_down <- ggplot(down,aes(x=down$Term,y=down$Count)) + geom_histogram(stat='identity',colour='black',fill='white') + geom_text(aes(label=down$Count),hjust=-0.12) +
  xlab("Gene Ontology") + ylab("Nombre de gènes") + ggtitle(" Gènes sous-exprimés") + 
  scale_fill_discrete(guide = guide_legend(reverse = T)) + coord_flip() + labs(fill = 'Definition')



mcols(res, use.names=TRUE)
overRepresented