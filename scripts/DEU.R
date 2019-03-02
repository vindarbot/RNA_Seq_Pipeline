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
library(pasilla)



dir = setwd("~/Desktop/Data")


pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)

head("/Library/Frameworks/R.framework/Versions/3.5/Resources/library/pasilla/extdata/treated3fb.txt")

## ----systemFileCheck,echo=FALSE,results='hide'--------------------------------
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

## ----loadDEXSeq---------------------------------------------------------------
inDir = system.file("extdata", package="pasilla")

countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)

countFiles = list.files("DEU/Counts", pattern=".tsv$", full.names=TRUE)

flattenedFile = list.files("Reference", pattern="reference.DEXSeq.gff", full.names=TRUE)

flattenedFile

## ----sampleTable--------------------------------------------------------------
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3", 
                 "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",  
                "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end", 
               "single-end", "single-end", "paired-end", "paired-end" ) )


sampleTable = data.frame(
  row.names = c("Col_1","Col_2","hon4_1","hon4_2","hon4_3"),
  condition = c("Col","Col","hon4","hon4","hon4"),
  libType = c("paired-end","paired-end","paired-end","paired-end","paired-end")
)


## ----makeecs, eval=TRUE-------------------------------------------------------
suppressPackageStartupMessages( library( "DEXSeq" ) )

dxd = DEXSeqDataSetFromHTSeq (
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

## ----start--------------------------------------------------------------------
genesForSubset = read.table( 
  file.path(inDir, "geneIDsinsubset.txt"), 
  stringsAsFactors=FALSE)[[1]]

dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

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
dxr1

## ----results2,cache=TRUE------------------------------------------------------
elementMetadata(dxr1)$description

## ----tallyExons---------------------------------------------------------------
table ( dxr1$padj < 0.1 )

## ----tallyGenes---------------------------------------------------------------
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )

## ----MvsA, dev='png', resolution=200------------------------------------------
plotMA( dxr1, cex=0.8 )

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
plotDEXSeq( dxr2, "AT1G01010", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

## ----checkClaim,echo=FALSE----------------------------------------------------
wh = (dxr2$groupID=="AT1G01010")
stopifnot(sum(dxr2$padj[wh] < formals(plotDEXSeq)$FDR)==1)

## ----plot2, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr2, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE,
            cex.axis=1.2, cex=1.3, lwd=2 )

## ----plot3, fig.height=8, fig.width=12----------------------------------------
plotDEXSeq( dxr2, "AT1G01010", expression=FALSE, norCounts=TRUE,
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

