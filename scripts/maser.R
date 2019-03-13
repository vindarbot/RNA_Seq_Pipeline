library(ggplot2)
library(tidyverse)
library(maser)
library(biomaRt)
### QC des données géénréas par rMATS

?maser()
dir = setwd("~/Desktop/Data")


path2 <- system.file("extdata", file.path("MATS_output"),
                    package = "maser")


path <- "rMATS_hon4HS/"

hon4 <- maser(path, c("ColHS","hon4_HS"), ftype="JCEC")

head(summary(hon4, type = "SE")[, 1:8])

hon4_filt <- filterByCoverage(hon4,avg_reads = 5)

hon4_top <- topEvents(hon4_filt, fdr = 0.05, deltaPSI = 0.1)

ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")
symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

merge <- merge(x=symbol,y=annotation(hon4_top),by.x="tair_locus",by.y="GeneID",sort=F,all.y=T)


table_hon4 <- annotation(hon4_top)


volcano(hon4,fdr=0.05,deltaPSI=0.1)

events <- geneEvents(hon4_filt,fdr=0.05,deltaPSI=0.1)

plotTranscripts(hon4,type="SE",gtf="Reference")

