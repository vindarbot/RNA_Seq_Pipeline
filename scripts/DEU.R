library(Rsamtools)
library(GenomicFeatures)
library(DESeq2)
library(GOstats)
library(Category)
library(org.At.tair.db)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(DEXSeq)
library(docopt)




dir = setwd("~/Desktop/Data")

pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )

gene_description <- read.delim("~/Desktop/RNA_SEQ/gene_description_20131231.txt",
                               header=FALSE, quote="")


## ----systemFileCheck,echo=FALSE,results='hide'--------------------------------
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

## ----loadDEXSeq---------------------------------------------------------------

countFiles = list.files("DEU/Counts", pattern=".tsv$", full.names=TRUE)

flattenedFile = list.files("Reference", pattern="reference.DEXSeq.gff", full.names=TRUE)


## ----sampleTable--------------------------------------------------------------

sampleTable = data.frame(
  row.names = c("Col_1","Col_2","hon4_1","hon4_2","hon4_3"),
  condition = c("Col","Col","hon4","hon4","hon4"),
  libType = c("paired-end","paired-end","paired-end","paired-end","paired-end")
)


## ----makeecs, eval=TRUE-------------------------------------------------------

dxd = DEXSeqDataSetFromHTSeq (
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )


## ----seeColData---------------------------------------------------------------
colData(dxd)

## ----seeCounts----------------------------------------------------------------
head( counts(dxd), 5 )

## ----seeSplitted--------------------------------------------------------------
split( seq_len(ncol(dxd)), colData(dxd)$exon )

## ----seeCounts2---------------------------------------------------------------
head( featureCounts(dxd), 5 )

## ----fData--------------------------------------------------------------------
head( rowData(dxd), 3 )

## ----pData--------------------------------------------------------------------
sampleAnnotation( dxd )

## ----sizeFactors1-------------------------------------------------------------
dxd = estimateSizeFactors( dxd )

## ----estDisp1-----------------------------------------------------------------
dxd = estimateDispersions( dxd )

## ----fitdiagnostics, dev='png', resolution=220--------------------------------
plotDispEsts( dxd )

## ----testForDEU1,cache=TRUE---------------------------------------------------
dxd = testForDEU( dxd )

## ----estFC,cache=TRUE---------------------------------------------------------
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

## ----results1,cache=TRUE------------------------------------------------------
dxr1 = DEXSeqResults( dxd )


genes_up <- dxr1[ which(dxr1$padj < 0.1 & dxr1$log2fold_hon4_Col > 0), ]
genes_up <- genes_up[order(genes_up$padj, decreasing = F),]

genes_down <- dxr1[ which(dxr1$padj < 0.1 & dxr1$log2fold_hon4_Col < 0), ]
genes_down <- genes_down[order(genes_up$padj, decreasing = F),]


genes_down <- as.data.frame(genes_down)
genes_up <- as.data.frame(genes_up)

genes_up <- genes_up[,c(1,6,7,10,15,16)]
genes_down <- genes_down[,c(1,6,7,10,15,16)]

identifiants_genes_up <- gsub(":[E][0-9]+$","",rownames(genes_up))
identifiants_genes_down <- gsub(":[E][0-9]+$","",rownames(genes_down))

genes_up$ID <- identifiants_genes_up
genes_down$ID <- identifiants_genes_down

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)

## Le fichier gene_description a été récupéré sur TAIR10 et nous permet d'accéder aux annotations
# des gènes, à partir de leur identifiants.


genenames = gsub("[.][1234567890]", "",
                 gene_description[,1])

gene_description[,1]=genenames

## Ici on récupère la description des gènes différentiellement exprimés
genes_match_rows <- match(genes_signif$groupID, gene_description[,1])

genes_match_up <- match((genes_up$groupID), gene_description[,1])

Code_to_description <- gene_description[genes_match_rows,c(1,3)]

Code_to_description_up <- gene_description[genes_match_up,c(1,3)]


colnames(Code_to_description) <- c("Code","Description")

colnames(Code_to_description_up) <- c("Code","Description")


genes_match_down <- match((genes_down$groupID), gene_description[,1])

Code_to_description_down <- gene_description[genes_match_down,c(1,3)]

colnames(Code_to_description_down) <- c("Code","Description")


genes_up$Description <- Code_to_description_up$Description

genes_down$Description <- Code_to_description_down$Description



ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")

TAIRID <- genes_down$ID

symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

genes_up <- (merge(x=symbol,y=genes_up,by.x="tair_locus",by.y="ID"))

genes_down <- (merge(x=symbol,y=genes_down,by.x="tair_locus",by.y="ID"))


write_tsv(as.data.frame(genes_up),"DEU/genes_up_DEXSEQ.txt")
write_tsv(as.data.frame(genes_down),"DEU/genes_down_DEXSEQ.txt")























## ----results2,cache=TRUE------------------------------------------------------
elementMetadata(dxr1)$description

## ----tallyExons---------------------------------------------------------------
table ( dxr1$padj < 0.1 )

## ----tallyGenes---------------------------------------------------------------
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )

## ----MvsA, dev='png', resolution=200------------------------------------------
plotMA( dxr1)

## ----design-------------------------------------------------------------------
sampleAnnotation(dxd)

## ----formulas2----------------------------------------------------------------
formulaFullModel    =  ~ sample + exon + libType:exon + condition:exon
formulaReducedModel =  ~ sample + exon + libType:exon 

## ----estDisps_again, cache=TRUE, results='hide'-------------------------------
dxd = estimateDispersions( dxd, formula = formulaFullModel )

## ----test_again, cache=TRUE---------------------------------------------------
dxd = testForDEU( dxd, 
                  reducedModel = formulaReducedModel, 
                  fullModel = formulaFullModel )

## ----res_again----------------------------------------------------------------
dxr2 = DEXSeqResults( dxd )

## ----table2-------------------------------------------------------------------
table( dxr2$padj < 0.1 )

## ----table3-------------------------------------------------------------------
table( before = dxr1$padj < 0.1, now = dxr2$padj < 0.1 )

## ----plot1, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr1, "AT1G64790", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

## ----checkClaim,echo=FALSE----------------------------------------------------
wh = (dxr2$groupID=="AT1G64790")
stopifnot(sum(dxr2$padj[wh] < formals(plotDEXSeq)$FDR)==1)

## ----plot2, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr1, "AT1G64790", displayTranscripts=TRUE, legend=TRUE,
            cex.axis=1.2, cex=1.3, lwd=2 )

## ----plot3, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr2, "AT1G64790", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

## ----plot4, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr2, "AT1G01010", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

## ----DEXSeqHTML,cache=TRUE, eval=FALSE----------------------------------------
DEXSeqHTML( dxr2, FDR=0.1, color=c("#FF000080", "#0000FF80") )

## ----para1,cache=TRUE,results='hide', eval=FALSE------------------------------
#  BPPARAM = MultiCoreParam(workers=4)
#  dxd = estimateSizeFactors( dxd )
#  dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
#  dxd = testForDEU( dxd, BPPARAM=BPPARAM)
#  dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)

## ----alldeu, cache=TRUE-------------------------------------------------------
dxr = DEXSeq(dxd)
class(dxr)

## ----buildExonCountSetLoadPacks,cache=TRUE, eval=FALSE------------------------
#  library(GenomicRanges)
#  library(GenomicFeatures)
#  library(GenomicAlignments)

## ----buildExonCountSetDownloadAnno,cache=TRUE, eval=FALSE---------------------
#  hse = makeTranscriptDbFromBiomart( biomart="ensembl",
#     dataset="hsapiens_gene_ensembl" )

## ----buildExonCountSetDisjoin,cache=TRUE, eval=FALSE--------------------------
#  exonicParts = disjointExons( hse, aggregateGenes=FALSE )

## ----buildExonCountSet2FindBAMs,cache=TRUE, eval=FALSE------------------------
#  bamDir = system.file( "extdata", package="parathyroidSE", mustWork=TRUE )
#  fls = list.files( bamDir, pattern="bam$", full=TRUE )

## ----buildExonCountSet2FindBAMs2,cache=TRUE, eval=FALSE-----------------------
#  bamlst = BamFileList( fls, index=character(), yieldSize=100000, obeyQname=TRUE )
#  SE = summarizeOverlaps( exonicParts, bamlst, mode="Union", singleEnd=FALSE,
#     ignore.strand=TRUE, inter.feature=FALSE, fragments=TRUE )

## ----buildExonCountSet3,cache=TRUE, eval=FALSE--------------------------------
#  colData(SE)$condition = c("A", "A", "B")
#  DEXSeqDataSetFromSE( SE, design= ~ sample + exon + condition:exon )

## ----acc----------------------------------------------------------------------
head( geneIDs(dxd) )
head( exonIDs(dxd) )

## ----grmethods----------------------------------------------------------------
interestingRegion = GRanges("chr2L", IRanges(start=3872658, end=3875302))
subsetByOverlaps( query=dxr, subject=interestingRegion )
findOverlaps( query=dxr, subject=interestingRegion )

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

