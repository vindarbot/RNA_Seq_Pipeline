library(ggplot2)
library(tidyverse)
library(maser)
library(biomaRt)
library(AnnotationHub)

### QC des données géénréas par rMATS


dir = setwd("~/Desktop/Data")


path <- "rMATS_hon4HS/"

hon4 <- maser(path, c("ColHS","hon4_HS"), ftype="JCEC")

head(summary(hon4, type = "SE")[, 1:8],50)

hon4_filt <- filterByCoverage(hon4,avg_reads = 5)

head(summary(hon4_filt, type = "SE")[, 1:8],50)



hon4_top <- topEvents(hon4_filt, fdr = 0.05, deltaPSI = 0.1)



ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")
symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

annotation(hon4_top) <- merge(x=symbol,y=annotation(hon4_top),by.x="tair_locus",by.y="GeneID",sort=F,all.y=T)

table_hon4 <- annotation(hon4_top)

gtf_path <- file.path("Reference/reference.gtf")
tair_gtf <- rtracklayer::import.gff(gtf_path)


hon4_mapped <- mapTranscriptsToEvents(hon4_top,tair_gtf )
annotation(hon4_mapped)

volcano(hon4,fdr=0.05,deltaPSI=0.1)

events <- geneEvents(hon4_filt,fdr=0.05,deltaPSI=0.1)

plotTranscripts(hon4,type="SE",gtf="Reference")


